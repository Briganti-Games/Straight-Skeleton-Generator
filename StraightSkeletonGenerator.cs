using Briganti.StraightSkeletons.Priority_Queue;
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Profiling;

namespace Briganti.StraightSkeletonGeneration
{
	public class StraightSkeletonGenerator : IEventTimeChangeListener
	{
		private WavefrontGraph wavefront;
		private StraightSkeleton straightSkeleton;

		private FastStructPriorityQueue<int> eventQueue;

		private int[] wavefrontToStraightSkeletonVertexIndices;

		private float time = 0.0f;
		private readonly float maxEventTime;

		private List<int> eventBatches = new List<int>();
		private int eventBatchIndex = 0;

		private List<int> edgeEventsAtSamePos = new List<int>();

		private bool debug;


		public StraightSkeletonGenerator(PolygonWithHoles polygonWithHoles, float maxEventTime, bool performDebugChecks)
		{
			this.debug = performDebugChecks;
			this.maxEventTime = maxEventTime;

			// create the initial graph based on the polygon
			wavefront = new WavefrontGraph(polygonWithHoles, this, debug);

			// generate the event queue
			eventQueue = new FastStructPriorityQueue<int>(wavefront.maxEdges);

			// we keep track of a mapping of vertex indices in the wavefront to vertex indices in the straight skeleton
			// this is because there might be multiple vertices in the same spot in the wavefront (for split events)
			// but we merge them in the straight skeleton.
			wavefrontToStraightSkeletonVertexIndices = new int[wavefront.maxVertices];
			for (int i = 0; i < wavefront.maxVertices; ++i)
			{
				wavefrontToStraightSkeletonVertexIndices[i] = -1; // not set!
			}


			// we copy all the initial vertices & edges to the straight skeleton - they are always part
			int nStraightSkeletonVertices = wavefront.nVertices + (wavefront.nVertices);

			// count the number of possible edges in the straight skeleton
			int nStartEdges = wavefront.nEdges; // each initial edge is one
			int nStraightSkeletonEdges = (2 * wavefront.nVertices) - 3; // then this is the max number that can occur because of edge events
			nStraightSkeletonEdges += wavefront.nEdges; // there is a worst-case chance that every edge will be reproduced again because no edge events occur before the maxEventTime
			nStraightSkeletonEdges *= 2; // we add edges in both directions to be able to generate polygon slabs from the straight skeleton more easily
			straightSkeleton = new StraightSkeleton(nStraightSkeletonVertices, nStraightSkeletonEdges + nStartEdges);
			for (int i = 0; i < wavefront.nVertices; ++i)
			{
				straightSkeleton.AddVertex(wavefront.vertices[i], 0);

				// the initial vertices map to the same position in the straight skeleton
				wavefrontToStraightSkeletonVertexIndices[i] = i;
			}
			for (int i = 0; i < wavefront.nEdges; ++i)
			{
				ref var edge = ref wavefront.edges[i];
				straightSkeleton.AddEdge(new Edge(edge.prevVertexIndex, edge.nextVertexIndex));
			}

			// initialize the events now - we are ready to receive them
			wavefront.GenerateInitialEvents();

		}

		public void Step()
		{
			if (!IsDone())
			{
				if (IsEventBatchesEmpty()) GenerateEventBatches();
				if (!IsEventBatchesEmpty()) ProcessEventBatch();
			}
		}

		public bool IsDone()
		{
			return eventQueue.Count == 0 && IsEventBatchesEmpty();
		}


		private bool IsEventBatchesEmpty()
		{
			return eventBatchIndex == eventBatches.Count;
		}

		public StraightSkeleton Generate()
		{
			// now we keep popping those events until we're done
			while (eventQueue.Count > 0)
			{
				GenerateEventBatches();
				while (!IsEventBatchesEmpty())
				{
					ProcessEventBatch();
				}
			}

			return straightSkeleton;
		}

		private void GenerateEventBatches()
		{
			if (eventBatchIndex != eventBatches.Count) throw new InvalidOperationException($"There are previous batches that still need to be resolved!");

			eventBatches.Clear();
			eventBatchIndex = 0;

			int edgeIndex = DequeueNextEvent();
			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];

			// we ran out of events in the queue - tap out early
			while (edgeEvent.eventType == EventType.None || edgeEvent.eventType == EventType.NotInWavefront)
			{
				if (eventQueue.Count == 0) return;
				edgeIndex = DequeueNextEvent();
				edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			}

			// if this is the first event that is over time, we just go over all remaining edges in the queue and spawn them at their position at maxEventTime
			if (edgeEvent.eventTime > maxEventTime + Geometry.EPS)
			{
				SpawnRemainingEdgesAtMaxTime(edgeIndex);
				return;
			}

			eventBatches.Add(edgeIndex);

			// find all events at the same time
			while (eventQueue.Count > 0 && IsAtSameTime(edgeIndex, eventQueue.First))
			{
				edgeIndex = DequeueNextEvent();
				edgeEvent = ref wavefront.edgeEvents[edgeIndex];
				if (edgeEvent.eventType == EventType.None) continue; // skip this event
				eventBatches.Add(edgeIndex);
			}

			// we now sort the events by their event pos - that way, we can create only one vertex for each set of edge events that end at the same position
			//eventBatches.Sort((v1, v2) => CompareEdgeEvent(v1, v2));
		}

		private int DequeueNextEvent()
		{
			int edgeIndex = eventQueue.Dequeue();
			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			edgeEvent.queueId = -1;

			while (edgeEvent.eventType == EventType.None && eventQueue.Count > 0)
			{
				edgeIndex = eventQueue.Dequeue();
				edgeEvent = ref wavefront.edgeEvents[edgeIndex];
				edgeEvent.queueId = -1;
			}

			// return the first proper event
			return edgeIndex;
		}

		private int CompareEdgeEvent(int edgeIndex1, int edgeIndex2)
		{
			ref var edgeEvent1 = ref wavefront.edgeEvents[edgeIndex1];
			ref var edgeEvent2 = ref wavefront.edgeEvents[edgeIndex2];

			int cmp = edgeEvent1.eventType.CompareTo(edgeEvent2.eventType);
			if (cmp == 0) cmp = Geometry.CompareTo(edgeEvent1.eventPos, edgeEvent2.eventPos);
			return cmp;
		}

		private void ProcessEventBatch()
		{

			// process a non-batch event
			ref var edgeEvent = ref wavefront.edgeEvents[eventBatches[eventBatchIndex]];
			if (!edgeEvent.eventType.IsBatchEvent())
			{
				ProcessNonBatchEvent(eventBatches[eventBatchIndex]);
				++eventBatchIndex;
			}

			// go over the batch and split them up in groups at the same pos
			else
			{
				edgeEventsAtSamePos.Clear();
				/*for (int i = eventBatchIndex; i < eventBatches.Count; ++i)
				{
					edgeEventsAtSamePos.Add(eventBatches[i]);
					if (i == eventBatches.Count - 1 || !IsSameTypeAtSamePos(eventBatches[i], eventBatches[i + 1]))
					{
						ProcessBatchEvents(edgeEventsAtSamePos);
						eventBatchIndex = i + 1;
						break;
					}
				}*/

				edgeEventsAtSamePos.Add(eventBatches[eventBatchIndex]);
				for (int i = eventBatchIndex + 1; i < eventBatches.Count; ++i)
				{
					if (IsSameTypeAtSamePos(eventBatches[eventBatchIndex], eventBatches[i]))
					{
						edgeEventsAtSamePos.Add(eventBatches[i]);
						eventBatches.RemoveAt(i);
						--i;
					}
				}
				ProcessBatchEvents(edgeEventsAtSamePos);
				++eventBatchIndex;
			}

			// if we emptied the entire set of events that happened at the same time, we update the wavefront
			time = edgeEvent.eventTime;
			wavefront.UpdateWavefront(time);

			// validate the graph to catch bugs early
			//Debug.Log("Wavefront consists of " + wavefront.nVertices + " vertices and " + wavefront.nEdges + " edges, while straight skeleton consists of " + straightSkeleton.nVertices + " vertices and " + straightSkeleton.nEdges + " edges");
			wavefront.ValidateState(time);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private bool IsAtSameTime(int edgeIndex1, int edgeIndex2)
		{
			ref EdgeEvent edgeEvent1 = ref wavefront.edgeEvents[edgeIndex1];
			ref EdgeEvent edgeEvent2 = ref wavefront.edgeEvents[edgeIndex2];
			return math.abs(edgeEvent1.eventTime - edgeEvent2.eventTime) < Geometry.EPS_LOWPRECISION;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private bool IsSameTypeAtSamePos(int edgeIndex1, int edgeIndex2)
		{
			ref EdgeEvent edgeEvent1 = ref wavefront.edgeEvents[edgeIndex1];
			ref EdgeEvent edgeEvent2 = ref wavefront.edgeEvents[edgeIndex2];
			if (edgeEvent1.eventType != edgeEvent2.eventType) return false;
			return Geometry.CompareToLowPrecision(edgeEvent1.eventPos, edgeEvent2.eventPos) == 0;
		}

		public void AddOrUpdateEdgeEvent(int edgeIndex)
		{
			Profiler.BeginSample("EventQueue.AddOrUpdateEdgeEvent");

			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];

			// could be removed and in the batch to be processed at this point
			if (eventQueue.Contains(edgeEvent.queueId))
			{
				// if this is a none-event, we either remove it or ignore it
				if (edgeEvent.eventType == EventType.None || edgeEvent.eventType == EventType.NotInWavefront)
				{
					eventQueue.Remove(edgeEvent.queueId);
					edgeEvent.queueId = -1;
				}

				// update the priority of the valid event
				else
				{
					eventQueue.UpdatePriority(edgeEvent.queueId, edgeEvent.eventTime);
				}
			}
			else
			{
				eventQueue.Enqueue(edgeIndex, edgeEvent.eventTime, out int queueId);
				edgeEvent.queueId = queueId;
			}

			Profiler.EndSample();
		}

		private void ProcessNonBatchEvent(int edgeIndex)
		{
			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			if (edgeEvent.eventType.IsBatchEvent()) throw new ArgumentException($"Event {edgeEvent} is a batch event that we are trying to process separately.");

			switch (edgeEvent.eventType)
			{
				case EventType.Nook:
					ProcessNookEvent(edgeIndex);
					break;
			}
		}

		private void ProcessBatchEvents(List<int> edgeIndices)
		{
			if (edgeIndices.Count == 0) throw new ArgumentException($"You always need at least one event to process.");

			int firstEdgeEventIndex = edgeIndices[0];
			ref EdgeEvent firstEdgeEvent = ref wavefront.edgeEvents[firstEdgeEventIndex];
			if (!firstEdgeEvent.eventType.IsBatchEvent()) throw new ArgumentException($"Event {firstEdgeEvent} is not a batch event but we are still trying to process {edgeIndices.Count} events at the same time!");

			// edge event - easy
			switch (firstEdgeEvent.eventType)
			{
				case EventType.Edge:
					ProcessEdgeEventsAtSamePos(edgeIndices);
					break;

				case EventType.Split:
					ProcessSplitEvent(edgeIndices);
					break;

				default:
					throw new NotSupportedException($"{firstEdgeEvent} is not a batch event.");
			}
		}

		public void ProcessEdgeEventsAtSamePos(List<int> eventEdgeIndices)
		{
			// we are going to be merging all vertices generated in the wavefront, because their difference don't matter to us
			int newSKVertexIndex = -1;
			float2 newVertex = float2.zero;
			for (int i = 0; i < eventEdgeIndices.Count; ++i)
			{
				int edgeIndex = eventEdgeIndices[i];

				ref var edgeEvent = ref wavefront.edgeEvents[edgeIndex];
				if (edgeEvent.eventType != EventType.Edge)
				{
					if (edgeEvent.eventType == EventType.Nook) ProcessNookEvent(edgeIndex);
					continue;
				}

				// if this is the last edge in this batch, we want to update the adjacent edge events and vertex velocities
				bool updateEdgeEvents = (i == eventEdgeIndices.Count - 1);

				// first add the new vertex to the wavefront and update the wavefront vertices
				int newVertexIndex = wavefront.AddVertexToWavefrontAndRemoveEdge(edgeIndex, updateEdgeEvents);

				// invalid state - we didn't remove an edge
				if (newVertexIndex == -1) continue;

				if (i == 0)
				{
					newVertex = wavefront.vertices[newVertexIndex];
					newSKVertexIndex = straightSkeleton.AddVertex(newVertex, wavefront.GetVertexTime(newVertexIndex));
				}

				// map the new vertex to the first one in this batch
				EnsureWavefrontMappingCapacity();
				wavefrontToStraightSkeletonVertexIndices[newVertexIndex] = newSKVertexIndex;

				// we now also add an edge to the straight skeleton
				ref var edge = ref wavefront.edges[edgeIndex];

				// see if this warrants an edge
				int prevSKVertexIndex = GetStraightSkeletonVertexAtCreationTime(edge.prevVertexIndex);
				int nextSKVertexIndex = GetStraightSkeletonVertexAtCreationTime(edge.nextVertexIndex);
				float2 prevSKVertex = straightSkeleton.vertices[prevSKVertexIndex];
				float2 nextSKVertex = straightSkeleton.vertices[nextSKVertexIndex];

				// add the edges to the straight skeleton
				if (math.distance(prevSKVertex, newVertex) > Geometry.EPS)
				{

					AddStraightSkeletonArc(prevSKVertexIndex, newSKVertexIndex);
				}
				if (math.distance(nextSKVertex, newVertex) > Geometry.EPS)
				{
					AddStraightSkeletonArc(newSKVertexIndex, nextSKVertexIndex);
				}
			}
		}

		private List<Edge> newStraightSkeletonEdges = new List<Edge>();
		private List<int> splitVertexIndices = new List<int>();
		private List<int> newVertexIndices = new List<int>();
		public void ProcessSplitEvent(List<int> eventEdgeIndices)
		{
			// we are going to be merging all vertices generated in the wavefront, because their difference don't matter to us
			newStraightSkeletonEdges.Clear();
			splitVertexIndices.Clear();
			newVertexIndices.Clear();

			// first, make sure this is still a valid spit - another split or event might have deleted some of the edges or involved vertices
			bool splitCancelled = false;
			for (int i = 0; i < eventEdgeIndices.Count; ++i)
			{
				int edgeIndex = eventEdgeIndices[i];
				if (!wavefront.IsPartOfWavefront(edgeIndex))
				{
					splitCancelled = true;
					break;
				}
				else
				{

					// it can also happen that the vertex that is causing the split, has been deleted by an edge event before - we ignore this split as well!
					ref var edgeEvent = ref wavefront.edgeEvents[edgeIndex];
					int vertexIndex = edgeEvent.reflexVertexIndex;
					if (vertexIndex == -1)
					{
						splitCancelled = true;
						break;
					}
				}
			}

			if (splitCancelled) return;

			// collect all vertices involved in the split, and generate one if the relevant edge is an edge split
			for (int i = 0; i < eventEdgeIndices.Count; ++i)
			{
				int edgeIndex = eventEdgeIndices[i];
				ref var edgeEvent = ref wavefront.edgeEvents[edgeIndex];
				int vertexIndex = edgeEvent.reflexVertexIndex;
				ref var vertexData = ref wavefront.vertexDatas[vertexIndex];

				// add the reflex vertex to the list
				if (vertexData.type != WavefrontVertexType.Reflex) throw new ArgumentException($"Split edge event {edgeEvent} was triggered by a non-reflex vertex {vertexData}.");
				splitVertexIndices.Add(vertexIndex);

				// now if we're splitting an edge, we first split the edge and add the split point to the list
				if (vertexData.splitPoint == SplitPoint.Edge)
				{
					int newVertexIndex = wavefront.SplitEdge(edgeIndex);
					splitVertexIndices.Add(newVertexIndex);
				}

				// if we're not splitting an edge for we have a split at one of the edge ends, we actually add the endpoint as well if it is NOT a relfex vertex
				// this is because when two reflex vertices collide, they will also register to each other as a split point, but will already have been added earlier
				else
				{
					ref var edge = ref wavefront.edges[edgeIndex];
					int splitVertexIndex = edge.GetPoint(vertexData.splitPoint);
					ref var splitVertexData = ref wavefront.vertexDatas[splitVertexIndex];
					if (splitVertexData.type == WavefrontVertexType.Convex && !splitVertexIndices.Contains(splitVertexIndex)) splitVertexIndices.Add(splitVertexIndex);
				}
			}

			// now that we have a bunch of points all at the same position - some reflex vertices, some just split vertices, some endpoints of an edge,
			// we can actually perform the split of the graph into parts.
			wavefront.SplitGraphAtVertices(splitVertexIndices, eventEdgeIndices[0], newStraightSkeletonEdges, newVertexIndices);

			// we map every newly created vertex at the split position to the same one in the straight skeleton
			EnsureWavefrontMappingCapacity();
			bool first = true;
			int newSKVertexIndex = -1;
			foreach (int vertexIndex in newVertexIndices)
			{
				if (first)
				{
					float2 newVertex = wavefront.vertices[vertexIndex];
					newSKVertexIndex = straightSkeleton.AddVertex(newVertex, wavefront.GetVertexTime(vertexIndex));
					first = false;
				}

				if (newSKVertexIndex == -1) throw new InvalidOperationException($"New SK Vertex should have been assigned by now.");
				wavefrontToStraightSkeletonVertexIndices[vertexIndex] = newSKVertexIndex;
			}

			// now add all the reflex arcs that were created to the skeleton
			AddReflexArcsToSkeleton(newStraightSkeletonEdges);
		}

		private int GetStraightSkeletonVertexAtCreationTime(int wavefrontVertexIndex)
		{
			ref VertexData vertexData = ref wavefront.vertexDatas[wavefrontVertexIndex];
			return GetStraightSkeletonVertexAtTime(wavefrontVertexIndex, vertexData.creationTime);
		}

		private int GetStraightSkeletonVertexAtMaxTime(int vertexIndex)
		{
			return GetStraightSkeletonVertexAtTime(vertexIndex, maxEventTime);
		}

		private int GetStraightSkeletonVertexAtTime(int vertexIndex, float time)
		{
			// if there is a mapping, this means we already forwarded the vertex,
			// OR it was spawned at the max event time - we return it as is
			if (wavefrontToStraightSkeletonVertexIndices[vertexIndex] != -1)
			{
				return wavefrontToStraightSkeletonVertexIndices[vertexIndex];
			}

			// otherwise, we forward it and spawn it
			else
			{
				ref VertexData vertexData = ref wavefront.vertexDatas[vertexIndex];

				float2 newVertex = wavefront.GetVertexPosAtTime(vertexIndex, time);
				int newSKVertexIndex = straightSkeleton.AddVertex(newVertex, time);
				wavefrontToStraightSkeletonVertexIndices[vertexIndex] = newSKVertexIndex;
				return newSKVertexIndex;
			}
		}

		private void AddReflexArcsToSkeleton(List<Edge> newEdges)
		{
			// nothing got split - we don't do anything, because the split got cancelled (probably because it occured in a previous simultaneous event)
			if (newEdges.Count == 0) return;

			// create one new vertex and map the rest to this new one
			EnsureWavefrontMappingCapacity();

			// for each reflex vertex, we draw an arc to the new split point
			for (int i = 0; i < newEdges.Count; ++i)
			{
				var edge = newEdges[i];

				int prevSKVertexIndex = GetStraightSkeletonVertexAtCreationTime(edge.prevVertexIndex);
				int nextSKVertexIndex = GetStraightSkeletonVertexAtCreationTime(edge.nextVertexIndex);
				if (prevSKVertexIndex == -1) throw new InvalidOperationException($"Wavefront has an invalid vertex {edge.prevVertexIndex} that wasn't added to the straight skeleton before!");
				if (nextSKVertexIndex == -1) throw new InvalidOperationException($"Wavefront has an invalid vertex {edge.nextVertexIndex} that wasn't added to the straight skeleton before!");

				// see if this warrants an edge
				float2 prevSKVertex = straightSkeleton.vertices[prevSKVertexIndex];
				float2 nextSKVertex = straightSkeleton.vertices[nextSKVertexIndex];

				// add the edges to the straight skeleton
				if (math.distance(prevSKVertex, nextSKVertex) > Geometry.EPS)
				{
					AddStraightSkeletonArc(prevSKVertexIndex, nextSKVertexIndex);

				}
			}
		}

		public void ProcessNookEvent(int eventEdgeIndex)
		{
			// get the event and do a sanity check
			ref var edgeEvent = ref wavefront.edgeEvents[eventEdgeIndex];
			ref var edge = ref wavefront.edges[eventEdgeIndex];
			if (edgeEvent.eventType != EventType.Nook) throw new InvalidOperationException($"There is no nook event associated with edge {edge}.");

			// remove the nook from the wavefront
			wavefront.RemoveNook(eventEdgeIndex);

			// this is easy - just spawn the edge, since the vertices already exist
			int prevSKVertexIndex = GetStraightSkeletonVertexAtTime(edge.prevVertexIndex, edgeEvent.eventTime);
			int nextSKVertexIndex = GetStraightSkeletonVertexAtTime(edge.nextVertexIndex, edgeEvent.eventTime);

			// see if this warrants an edge
			float2 prevSKVertex = straightSkeleton.vertices[prevSKVertexIndex];
			float2 nextSKVertex = straightSkeleton.vertices[nextSKVertexIndex];

			// add the edges to the straight skeleton
			if (math.distance(prevSKVertex, nextSKVertex) > Geometry.EPS)
			{
				AddStraightSkeletonArc(prevSKVertexIndex, nextSKVertexIndex);
			}
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private void AddStraightSkeletonArc(int skVertexIndex1, int skVertexIndex2)
		{
			if (skVertexIndex1 == -1) throw new ArgumentException($"Trying to generate an arc from non-existing vertex {skVertexIndex1}");
			if (skVertexIndex2 == -1) throw new ArgumentException($"Trying to generate an arc from non-existing vertex {skVertexIndex2}");
			straightSkeleton.AddEdge(skVertexIndex1, skVertexIndex2);
			straightSkeleton.AddEdge(skVertexIndex2, skVertexIndex1);

		}

		private void SpawnRemainingEdgesAtMaxTime(int firstOverTimeEdgeIndex)
		{
			EnsureWavefrontMappingCapacity();

			// we generate straight skeleton edges from the current vertex to its position at the max time,
			// and also update the mapping from the old SK vertex position to the new one
			for (int i = 0; i < wavefront.nVertices; ++i)
			{
				ref VertexData vertexData = ref wavefront.vertexDatas[i];
				if (!vertexData.inWavefront) continue;
				SpawnEdgeFromPosToMaxTimePos(i);
			}

			// re-use this list for convenience
			SpawnEdgeAtMaxTime(firstOverTimeEdgeIndex);
			while (eventQueue.Count > 0)
			{
				var edgeIndex = eventQueue.Dequeue();
				SpawnEdgeAtMaxTime(edgeIndex);
			}
		}

		private void SpawnEdgeAtMaxTime(int edgeIndex)
		{
			ref Edge edge = ref wavefront.edges[edgeIndex];

			// if the vertex was spawned in the past, we spawn an up-to-date version at the max time
			int prevSKVertexIndexAtMaxTime = GetStraightSkeletonVertexAtMaxTime(edge.prevVertexIndex);
			int nextSKVertexIndexAtMaxTime = GetStraightSkeletonVertexAtMaxTime(edge.nextVertexIndex);
			AddStraightSkeletonArc(prevSKVertexIndexAtMaxTime, nextSKVertexIndexAtMaxTime);

		}

		private void SpawnEdgeFromPosToMaxTimePos(int vertexIndex)
		{
			ref VertexData vertexData = ref wavefront.vertexDatas[vertexIndex];
			if (vertexData.creationTime < maxEventTime - Geometry.EPS)
			{
				int startSKVertexIndex = GetStraightSkeletonVertexAtCreationTime(vertexIndex);
				wavefrontToStraightSkeletonVertexIndices[vertexIndex] = -1;
				int endSKVertexIndex = GetStraightSkeletonVertexAtMaxTime(vertexIndex);
				AddStraightSkeletonArc(startSKVertexIndex, endSKVertexIndex);
			}
		}

		public void DrawGizmos()
		{
			Gizmos.color = Color.red;

			for (int i = 0; i < straightSkeleton.nEdges; ++i)
			{
				var edge = straightSkeleton.edges[i];
				var prevVertex = straightSkeleton.vertices[edge.prevVertexIndex];
				var nextVertex = straightSkeleton.vertices[edge.nextVertexIndex];

				Gizmos.DrawLine(new Vector3(prevVertex.x, straightSkeleton.vertexTimes[edge.prevVertexIndex], prevVertex.y), new Vector3(nextVertex.x, straightSkeleton.vertexTimes[edge.nextVertexIndex], nextVertex.y));
			}

			if (!IsDone()) wavefront.DrawGizmos(time);
		}

		private void EnsureWavefrontMappingCapacity()
		{
			if (wavefrontToStraightSkeletonVertexIndices.Length != wavefront.maxVertices)
			{
				int oldLength = wavefrontToStraightSkeletonVertexIndices.Length;
				int[] newWavefrontToStraightSkeletonVertexIndices = new int[wavefront.maxVertices];
				Array.Copy(wavefrontToStraightSkeletonVertexIndices, newWavefrontToStraightSkeletonVertexIndices, oldLength);
				wavefrontToStraightSkeletonVertexIndices = newWavefrontToStraightSkeletonVertexIndices;

				for (int i = oldLength; i < wavefront.maxVertices; ++i)
				{
					wavefrontToStraightSkeletonVertexIndices[i] = -1;
				}
			}
		}
	}
}