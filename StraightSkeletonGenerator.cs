using Briganti.StraightSkeletons.Priority_Queue;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Profiling;
using UnityEngine.Rendering;
using UnityEngine.UI;

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


		public StraightSkeletonGenerator(PolygonWithHoles polygonWithHoles, float maxEventTime)
		{
			this.maxEventTime = maxEventTime;

			// create the initial graph based on the polygon
			wavefront = new WavefrontGraph(polygonWithHoles, this);

			// generate the event queue
			eventQueue = new FastStructPriorityQueue<int>(wavefront.maxVertices);

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
			int nStraightSkeletonEdges = wavefront.nEdges; // each initial edge is one
			nStraightSkeletonEdges += (2 * wavefront.nVertices) - 3; // then this is the max number that can occur because of edge events
			nStraightSkeletonEdges += wavefront.nEdges; // there is a worst-case chance that every edge will be reproduced again because no edge events occur before the maxEventTime
			straightSkeleton = new StraightSkeleton(nStraightSkeletonVertices, nStraightSkeletonEdges);
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
			if (eventQueue.Count > 0)
			{
				if (IsEventBatchesEmpty()) GenerateEventBatches();
				ProcessEventBatch();
			}
		}

		public bool IsDone()
		{
			return eventQueue.Count == 0;
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

			int edgeIndex = eventQueue.Dequeue();
			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			if (edgeEvent.eventType == EventType.None) return; // skip this event

			// if this is the first event that is over time, we just go over all remaining edges in the queue and spawn them at their position at maxEventTime
			if (edgeEvent.eventTime > maxEventTime)
			{
				SpawnRemainingEdgesAtMaxTime(edgeIndex);
				return;
			}

			time = edgeEvent.eventTime;
			eventBatches.Add(edgeIndex);

			// find all events at the same time
			while (eventQueue.Count > 0 && IsAtSameTime(edgeIndex, eventQueue.First))
			{
				edgeIndex = eventQueue.Dequeue();
				edgeEvent = ref wavefront.edgeEvents[edgeIndex];
				if (edgeEvent.eventType == EventType.None) continue; // skip this event
				eventBatches.Add(edgeIndex);
			}

			// at this point, we either gathered all batch events that occur at the same time, OR we bumped upon a non-batch event, and we need to process that one
			/*if (!edgeEvent.eventType.IsBatchEvent())
			{
				for (int i = 0; i < eventBatches.Count - 1; ++i)
				{
					edgeIndex = eventBatches[i];
					edgeEvent = ref wavefront.edgeEvents[edgeIndex];
					eventQueue.Enqueue(edgeIndex, edgeEvent.eventTime, out int queueId);
					edgeEvent.queueId = queueId;
				}

				eventBatches.Clear();
				eventBatches.Add(edgeIndex);
			}*/

			// we now sort the events by their event pos - that way, we can create only one vertex for each set of edge events that end at the same position
			eventBatches.Sort((v1, v2) => Geometry.CompareTo(v1, v2));

			// we start at the start
			eventBatchIndex = 0;
		}

		private void ProcessEventBatch()
		{
			// go over the batch and split them up in groups at the same pos
			edgeEventsAtSamePos.Clear();
			for (int i = eventBatchIndex; i < eventBatches.Count; ++i)
			{
				edgeEventsAtSamePos.Add(eventBatches[i]);

				if (i == eventBatches.Count - 1 || !IsAtSamePos(eventBatches[i], eventBatches[i + 1]))
				{
					ProcessEvents(edgeEventsAtSamePos);
					eventBatchIndex = i + 1;
					break;
				}
			}
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private bool IsAtSameTime(int edgeIndex1, int edgeIndex2)
		{
			ref EdgeEvent edgeEvent1 = ref wavefront.edgeEvents[edgeIndex1];
			ref EdgeEvent edgeEvent2 = ref wavefront.edgeEvents[edgeIndex2];
			return math.abs(edgeEvent1.eventTime - edgeEvent2.eventTime) < Geometry.EPS;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private bool IsAtSamePos(int edgeIndex1, int edgeIndex2)
		{
			ref EdgeEvent edgeEvent1 = ref wavefront.edgeEvents[edgeIndex1];
			ref EdgeEvent edgeEvent2 = ref wavefront.edgeEvents[edgeIndex2];
			return Geometry.CompareTo(edgeEvent1.eventPos, edgeEvent2.eventPos) == 0;
		}

		private bool IsEdgeEvent(int edgeIndex)
		{
			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			return edgeEvent.eventType == EventType.Edge;
		}

		public void AddEdgeEvent(int edgeIndex)
		{
			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			eventQueue.Enqueue(edgeIndex, edgeEvent.eventTime, out int queueId);
			edgeEvent.queueId = queueId;
		}

		public void UpdateEdgeEvent(int edgeIndex)
		{
			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];

			// could be removed and in the batch to be processed at this point
			if (eventQueue.Contains(edgeEvent.queueId))
			{
				eventQueue.UpdatePriority(edgeEvent.queueId, edgeEvent.eventTime);
			}
		}

		private void ProcessEvents(List<int> edgeEventIndices)
		{
			if (edgeEventIndices.Count == 0) throw new ArgumentException($"You always need at least one event to process.");

			int firstEdgeEventIndex = edgeEventIndices[0];
			ref EdgeEvent firstEdgeEvent = ref wavefront.edgeEvents[firstEdgeEventIndex];

			//if (!firstEdgeEvent.eventType.IsBatchEvent() && edgeEventIndices.Count > 1) throw new ArgumentException($"Event {firstEdgeEvent} is not a batch event but we are still trying to process {edgeEventIndices.Count} events at the same time!");

			// edge event - easy
			switch (firstEdgeEvent.eventType)
			{
				case EventType.Edge:
					ProcessEdgeEventsAtSamePos(edgeEventIndices);
					break;

				case EventType.Nook:
					ProcessNookEvent(firstEdgeEventIndex);
					break;
			}
		}

		public void ProcessEdgeEventsAtSamePos(List<int> eventEdgeIndices)
		{
			// we are going to be merging all vertices generated in the wavefront, because their difference don't matter to us
			int newSKVertexIndex = 0;
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
				if (i == 0)
				{
					newVertex = wavefront.vertices[newVertexIndex];
					newSKVertexIndex = straightSkeleton.AddVertex(newVertex, wavefront.GetVertexTime(newVertexIndex));
				}

				// spawn the vertex and store the mapping
				wavefrontToStraightSkeletonVertexIndices[newVertexIndex] = newSKVertexIndex;

				// we now also add two edges to the straight skeleton, from the old vertices to the new one, but we go through the mapping first
				ref var edge = ref wavefront.edges[edgeIndex];

				int prevSKVertexIndex = wavefrontToStraightSkeletonVertexIndices[edge.prevVertexIndex];
				int nextSKVertexIndex = wavefrontToStraightSkeletonVertexIndices[edge.nextVertexIndex];

				// see if this warrants an edge
				float2 prevSKVertex = straightSkeleton.vertices[prevSKVertexIndex];
				float2 nextSKVertex = straightSkeleton.vertices[nextSKVertexIndex];

				// add the edges to the straight skeleton
				if (math.distance(prevSKVertex, newVertex) > Geometry.EPS) straightSkeleton.AddEdge(prevSKVertexIndex, newSKVertexIndex);
				if (math.distance(nextSKVertex, newVertex) > Geometry.EPS) straightSkeleton.AddEdge(newSKVertexIndex, nextSKVertexIndex);
			}
		}

		public void ProcessNookEvent(int eventEdgeIndex)
		{
			ref var edgeEvent = ref wavefront.edgeEvents[eventEdgeIndex];
			ref var edge = ref wavefront.edges[eventEdgeIndex];

			if (edgeEvent.eventType != EventType.Nook) throw new InvalidOperationException($"There is no nook event associated with edge {edge}.");

			// this is easy - just spawn the edge, since the vertices already exist
			int prevSKVertexIndex = wavefrontToStraightSkeletonVertexIndices[edge.prevVertexIndex];
			int nextSKVertexIndex = wavefrontToStraightSkeletonVertexIndices[edge.nextVertexIndex];

			straightSkeleton.AddEdge(prevSKVertexIndex, nextSKVertexIndex);
		}

		private void SpawnRemainingEdgesAtMaxTime(int firstOverTimeEdgeIndex)
		{
			// we generate straight skeleton edges from the current vertex to its position at the max time,
			// and also update the mapping from the old SK vertex position to the new one
			for (int i = 0; i < wavefront.nVertices; ++i)
			{
				ref VertexData vertexData = ref wavefront.vertexDatas[i];
				if (!vertexData.inWavefront) continue;
				SpawnEdgeFromPosToMaxTimePos(i);
			}

			// we reset the wavefront mapping for each existing vertex that wasn't spawned at the exact max event time
			// because we want to generate new "forwarded" vertices for them!
			/*for (int i = 0; i < wavefront.nVertices; ++i)
			{
				ref VertexData vertexData = ref wavefront.vertexDatas[i];
				if (vertexData.creationTime < maxEventTime - Geometry.EPS)
				{
					wavefrontToStraightSkeletonVertexIndices[i] = -1;
				}
			}*/

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
			straightSkeleton.AddEdge(prevSKVertexIndexAtMaxTime, nextSKVertexIndexAtMaxTime);

		}

		private void SpawnEdgeFromPosToMaxTimePos(int vertexIndex)
		{
			ref VertexData vertexData = ref wavefront.vertexDatas[vertexIndex];
			if (vertexData.creationTime < maxEventTime - Geometry.EPS)
			{
				int startSKVertexIndex = wavefrontToStraightSkeletonVertexIndices[vertexIndex];
				wavefrontToStraightSkeletonVertexIndices[vertexIndex] = -1;
				int endSKVertexIndex = GetStraightSkeletonVertexAtMaxTime(vertexIndex);
				straightSkeleton.AddEdge(startSKVertexIndex, endSKVertexIndex);
			}
		}

		private int GetStraightSkeletonVertexAtMaxTime(int vertexIndex)
		{
			ref VertexData vertexData = ref wavefront.vertexDatas[vertexIndex];

			// if there is a mapping, this means we already forwarded the vertex,
			// OR it was spawned at the max event time - we return it as is
			if (wavefrontToStraightSkeletonVertexIndices[vertexIndex] != -1)
			{
				return wavefrontToStraightSkeletonVertexIndices[vertexIndex];
			}

			// otherwise, we forward it and spawn it
			else
			{
				float2 newVertex = wavefront.GetVertexPosAtTime(vertexIndex, maxEventTime);
				int newSKVertexIndex = straightSkeleton.AddVertex(newVertex, maxEventTime);
				wavefrontToStraightSkeletonVertexIndices[vertexIndex] = newSKVertexIndex;
				return newSKVertexIndex;
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

	}
}