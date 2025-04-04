using Briganti.StraightSkeletons.Priority_Queue;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Profiling;

namespace Briganti.StraightSkeletonGeneration
{
	public class StraightSkeletonGenerator : IEventTimeChangeListener
	{
		private struct QueueEvent
		{
			public int index;
			public EventType eventType;
			public float eventTime;
			public float2 eventPos;

			public int vertexIndex => QueueIndexToVertexIndex(index);
			public int edgeIndex => index;

			public QueueEvent(int edgeIndex, in EdgeEvent edgeEvent)
			{
				index = edgeIndex;
				eventType = edgeEvent.eventType;
				eventTime = edgeEvent.eventTime;
				eventPos = edgeEvent.eventPos;
			}

			public QueueEvent(int vertexIndex, in VertexData vertexData)
			{
				index = vertexIndex;
				eventType = EventType.VertexSplit;
				eventTime = vertexData.splitTime;
				eventPos = vertexData.splitPos;
			}

			[MethodImpl(MethodImplOptions.AggressiveInlining)]
			public static int VertexIndexToQueueIndex(int vertexIndex)
			{
				return -(vertexIndex + 1);
			}

			[MethodImpl(MethodImplOptions.AggressiveInlining)]
			public static int QueueIndexToVertexIndex(int queueIndex)
			{
				return -queueIndex - 1;
			}

			[MethodImpl(MethodImplOptions.AggressiveInlining)]
			public static int EdgeIndexToQueueIndex(int edgeIndex)
			{
				return edgeIndex;
			}

			public static readonly QueueEvent Invalid = new() { eventType = EventType.None };

			public override string ToString()
			{
				if (eventType == EventType.None) return "No event";
				else if (eventType == EventType.VertexSplit) return $"{eventType} event at time {eventTime}, pos {eventPos}, vertex index {vertexIndex}";
				else return $"{eventType} event at time {eventTime}, pos {eventPos}, edge index {edgeIndex}";
			}
		}

		private WavefrontGraph wavefront;
		private StraightSkeleton straightSkeleton;

		private FastStructPriorityQueue<int> eventQueue;

		private int[] wavefrontToStraightSkeletonVertexIndices;

		private float time = 0.0f;
		private readonly float maxEventTime;

		private List<QueueEvent> eventBatches = new List<QueueEvent>();
		private int eventBatchIndex = 0;

		private List<QueueEvent> queueEventsAtTheSamePos = new List<QueueEvent>();

		private bool debug;
		private bool verbose = false;


		public StraightSkeletonGenerator(PolygonWithHoles polygonWithHoles, float maxEventTime, bool performDebugChecks)
		{
			this.debug = performDebugChecks;
			this.maxEventTime = maxEventTime;

			// create the initial graph based on the polygon
			wavefront = new WavefrontGraph(polygonWithHoles, this, debug);

			// generate the event queue
			eventQueue = new FastStructPriorityQueue<int>(wavefront.maxEdges + wavefront.maxVertices);

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

		public void SetVerbose(bool verbose)
		{
			this.verbose = verbose;
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
			// events can be removed from the event batch, so we can suddenly be further than the number of eventbatches
			return eventBatchIndex >= eventBatches.Count;
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

			QueueEvent queueEvent = DequeueNextEvent();
			while (queueEvent.eventType == EventType.None && eventQueue.Count > 0) queueEvent = DequeueNextEvent();
			if (queueEvent.eventType == EventType.None) return;

			QueueEvent firstQueueEvent = queueEvent;

			// if this is the first event that is over time, we just go over all remaining edges in the queue and spawn them at their position at maxEventTime
			if (queueEvent.eventTime > maxEventTime + Geometry.EPS)
			{
				SpawnRemainingEdgesAtMaxTime(queueEvent);
				return;
			}

			eventBatches.Add(queueEvent);

			// find all events at the same time
			while (eventQueue.Count > 0 && IsAtSameTime(firstQueueEvent, eventQueue.First))
			{
				queueEvent = DequeueNextEvent();
				if (queueEvent.eventType == EventType.None) continue; // skip this event
				eventBatches.Add(queueEvent);
			}
		}

		private QueueEvent DequeueNextEvent()
		{
			int index = eventQueue.Dequeue();
			if (index >= 0)
			{
				int edgeIndex = index;
				ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
				if (!edgeEvent.eventType.IsValid()) return QueueEvent.Invalid;
				edgeEvent.queueId = -1;
				return new QueueEvent(index, edgeEvent);
			}
			else
			{
				int vertexIndex = QueueEvent.QueueIndexToVertexIndex(index);
				ref VertexData vertexData = ref wavefront.vertexDatas[vertexIndex];
				if (!vertexData.partOfSplitEvent) return QueueEvent.Invalid;
				vertexData.queueId = -1;
				return new QueueEvent(index, vertexData);
			}
		}

		private void ProcessEventBatch()
		{

			// process a non-batch event
			var queueEvent = eventBatches[eventBatchIndex];
			if (!queueEvent.eventType.IsBatchEvent())
			{
				++eventBatchIndex;
				ProcessNonBatchEvent(queueEvent);
			}

			// go over the batch and split them up in groups at the same pos
			else
			{
				queueEventsAtTheSamePos.Clear();
				queueEventsAtTheSamePos.Add(eventBatches[eventBatchIndex]);
				for (int i = eventBatchIndex + 1; i < eventBatches.Count; ++i)
				{
					if (IsSameTypeAtSamePos(eventBatches[eventBatchIndex], eventBatches[i]))
					{
						queueEventsAtTheSamePos.Add(eventBatches[i]);
						eventBatches.RemoveAt(i);
						--i;
					}
				}
				++eventBatchIndex;
				ProcessBatchEvents(queueEventsAtTheSamePos);
			}

			// if we emptied the entire set of events that happened at the same time, we update the wavefront
			time = queueEvent.eventTime;
			wavefront.UpdateWavefront(time);

			// validate the graph to catch bugs early
			//Debug.Log("Wavefront consists of " + wavefront.nVertices + " vertices and " + wavefront.nEdges + " edges, while straight skeleton consists of " + straightSkeleton.nVertices + " vertices and " + straightSkeleton.nEdges + " edges");
			wavefront.ValidateState(time);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private bool IsAtSameTime(in QueueEvent queueEvent, int queueIndex)
		{
			return math.abs(queueEvent.eventTime - GetEventTime(queueIndex)) < Geometry.EPS_LOWPRECISION;
		}

		private float GetEventTime(int queueIndex)
		{
			if (queueIndex >= 0) return wavefront.edgeEvents[queueIndex].eventTime;
			else return wavefront.vertexDatas[-(queueIndex + 1)].splitTime;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private bool IsSameTypeAtSamePos(QueueEvent queueEvent1, QueueEvent queueEvent2)
		{
			if (queueEvent1.eventType != queueEvent2.eventType) return false;
			return Geometry.CompareToLowPrecision(queueEvent1.eventPos, queueEvent2.eventPos) == 0;
		}

		public void AddOrUpdateEdgeEvent(int edgeIndex)
		{
			AddOrUpdateEdgeEvent(edgeIndex, true);
		}

		public void AddOrUpdateVertexEvent(int vertexIndex)
		{
			AddOrUpdateVertexEvent(vertexIndex, true);
		}

		public void AddOrUpdateEdgeEvent(int edgeIndex, bool validateEventBatches = true)
		{

			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			AddOrUpdateEvent(QueueEvent.EdgeIndexToQueueIndex(edgeIndex), edgeEvent.eventTime, edgeEvent.queueId, edgeEvent.eventType.IsValid(), validateEventBatches, out int newQueueId);
			edgeEvent.queueId = newQueueId;
		}

		public void AddOrUpdateVertexEvent(int vertexIndex, bool validateEventBatches = true)
		{

			ref VertexData vertexData = ref wavefront.vertexDatas[vertexIndex];
			bool valid = (vertexData.inWavefront && vertexData.partOfSplitEvent);
			AddOrUpdateEvent(QueueEvent.VertexIndexToQueueIndex(vertexIndex), vertexData.splitTime, vertexData.queueId, valid, validateEventBatches, out int newQueueId);
			vertexData.queueId = newQueueId;
		}

		private void AddOrUpdateEvent(int eventIndex, float eventTime, int queueId, bool validEvent, bool validateEventBatches, out int newQueueId)
		{
			Profiler.BeginSample("EventQueue.AddOrUpdateEvent.ValidateEventBatches");

			// if we are in the current event batch and we moved back in time from the current time, we remove ourselves from the batch
			// this means that a previous event in the batch has changed our event time, and we are not ready to be handled just yet.
			// this can ONLY happen if we are not in a queue, because everything in eventBatches was removed from the queue.
			if (queueId == -1)
			{
				int eventBatchIndex = eventBatches.FindIndex(queueEvent => queueEvent.index == eventIndex);
				if (eventBatchIndex != -1 && eventBatchIndex >= this.eventBatchIndex && eventTime >= time + Geometry.EPS_LOWPRECISION)
				{
					eventBatches.RemoveAt(eventBatchIndex);
				}

				// if a new or updated event should be in the batch right now but isn't, we totally invalidate the entire batch - we need to start all over again!
				if (eventBatches.Count > 0 && eventBatchIndex == -1 && IsAtSameTime(eventBatches[0], eventIndex))
				{
					Debug.Log("Clear " + string.Join(", ", eventBatches.Skip(this.eventBatchIndex)) + " because the new event " + eventIndex + " at " + eventTime + " should be part of this batch");
					for (int i = this.eventBatchIndex; i < eventBatches.Count; ++i)
					{
						QueueEvent queueEvent = eventBatches[i];

						// this is the actual same event that we're adding right now - don't add it yet, as it will be added later
						if (queueEvent.index == eventIndex) continue;

						if (queueEvent.eventType != EventType.None && queueEvent.eventType != EventType.NotInWavefront)
						{
							if (queueEvent.eventType == EventType.VertexSplit) AddOrUpdateVertexEvent(queueEvent.vertexIndex, false);
							else AddOrUpdateEdgeEvent(queueEvent.edgeIndex, false);
						}
					}
					eventBatches.Clear();
					this.eventBatchIndex = 0;
				}
			}

			Profiler.EndSample();

			Profiler.BeginSample("EventQueue.AddOrUpdateEdgeEvent");

			// could be removed and in the batch to be processed at this point
			if (eventQueue.Contains(queueId))
			{
				// if this is a none-event, we either remove it or ignore it
				if (!validEvent)
				{
					eventQueue.Remove(queueId);
					newQueueId = -1;
				}

				// update the priority of the valid event
				else
				{
					eventQueue.UpdatePriority(queueId, eventTime);
					newQueueId = queueId;
				}
			}
			else
			{
				if (validEvent)
				{
					eventQueue.Enqueue(eventIndex, eventTime, out newQueueId);
				}
				else
				{
					newQueueId = -1;
				}
			}

			Profiler.EndSample();
		}

		private void ProcessNonBatchEvent(QueueEvent queueEvent)
		{
			//ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			if (queueEvent.eventType.IsBatchEvent()) throw new ArgumentException($"Event {queueEvent} is a batch event that we are trying to process separately.");

			switch (queueEvent.eventType)
			{
				case EventType.Nook:
					ProcessNookEvent(queueEvent.edgeIndex);
					break;
			}
		}

		private void ProcessBatchEvents(List<QueueEvent> queueEvents)
		{
			if (queueEvents.Count == 0) throw new ArgumentException($"You always need at least one event to process.");

			QueueEvent firstQueueEvent = queueEvents[0];
			if (!firstQueueEvent.eventType.IsBatchEvent()) throw new ArgumentException($"Event {firstQueueEvent} is not a batch event but we are still trying to process {queueEvents.Count} events at the same time!");

			// edge event - easy
			switch (firstQueueEvent.eventType)
			{
				case EventType.Edge:
					ProcessEdgeEventsAtSamePos(queueEvents);
					break;

				case EventType.VertexSplit:
					ProcessSplitEvent(queueEvents);
					break;

				default:
					throw new NotSupportedException($"{firstQueueEvent} is not a batch event.");
			}
		}

		private void ProcessEdgeEventsAtSamePos(List<QueueEvent> queueEvents)
		{
			// we are going to be merging all vertices generated in the wavefront, because their difference don't matter to us
			int newSKVertexIndex = -1;
			float2 newVertex = float2.zero;
			for (int i = 0; i < queueEvents.Count; ++i)
			{
				int edgeIndex = queueEvents[i].edgeIndex;

				ref var edgeEvent = ref wavefront.edgeEvents[edgeIndex];
				if (edgeEvent.eventType != EventType.Edge)
				{
					if (edgeEvent.eventType == EventType.Nook) ProcessNookEvent(edgeIndex);
					continue;
				}

				// if this is the last edge in this batch, we want to update the adjacent edge events and vertex velocities
				bool updateEdgeEvents = (i == queueEvents.Count - 1);

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
		private void ProcessSplitEvent(List<QueueEvent> queueEvents)
		{
			// we are going to be merging all vertices generated in the wavefront, because their difference don't matter to us
			newStraightSkeletonEdges.Clear();
			splitVertexIndices.Clear();
			newVertexIndices.Clear();

			// first, make sure this is still a valid spit - another split or event might have deleted some of the edges or involved vertices
			bool splitCancelled = false;
			for (int i = 0; i < queueEvents.Count; ++i)
			{
				int vertexIndex = queueEvents[i].vertexIndex;
				if (!wavefront.vertexDatas[vertexIndex].partOfSplitEvent)
				{
					splitCancelled = true;
					break;
				}
			}

			if (splitCancelled) return;

			// collect all vertices involved in the split, and generate one if the relevant edge is an edge split
			float splitTime = 0;
			float2 splitPos = float2.zero;
			for (int i = 0; i < queueEvents.Count; ++i)
			{
				int vertexIndex = queueEvents[i].vertexIndex;
				ref var vertexData = ref wavefront.vertexDatas[vertexIndex];
				int edgeIndex = vertexData.splitEdge;

				// they are the same for all events because they're in the same batch
				if (i == 0)
				{
					splitTime = vertexData.splitTime;
					splitPos = vertexData.splitPos;
				}

				// add the reflex vertex to the list
				if (vertexData.type != WavefrontVertexType.Reflex) throw new ArgumentException($"Split edge event {vertexData} was triggered by a non-reflex vertex.");
				splitVertexIndices.Add(vertexIndex);

				// now if we're splitting an edge, we first split the edge and add the split point to the list
				// ... but never split the same edge twice!
				if (vertexData.splitPoint == SplitPoint.Edge)
				{
					if (wavefront.edgeEvents[edgeIndex].eventType != EventType.NotInWavefront)
					{
						int newVertexIndex = wavefront.SplitEdge(vertexIndex, edgeIndex, splitTime, splitPos);
						splitVertexIndices.Add(newVertexIndex);
					}
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
			wavefront.SplitGraphAtVertices(splitVertexIndices, splitTime, splitPos, newStraightSkeletonEdges, newVertexIndices);

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

		private void SpawnRemainingEdgesAtMaxTime(QueueEvent queueEvent)
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
			if (queueEvent.index >= 0) SpawnEdgeAtMaxTime(queueEvent.edgeIndex);
			while (eventQueue.Count > 0)
			{
				var queueIndex = eventQueue.Dequeue();
				if (queueIndex >= 0) SpawnEdgeAtMaxTime(queueIndex);
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