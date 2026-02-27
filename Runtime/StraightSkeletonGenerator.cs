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
				index = VertexIndexToQueueIndex(vertexIndex);
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

			[MethodImpl(MethodImplOptions.AggressiveInlining)]
			public bool IsEdgeEvent()
			{
				return index >= 0;
			}

			[MethodImpl(MethodImplOptions.AggressiveInlining)]
			public bool IsVertexEvent()
			{
				return !IsEdgeEvent();
			}

			public static readonly QueueEvent Invalid = new() { eventType = EventType.None };

			public override string ToString()
			{
				if (IsEdgeEvent()) return $"{eventType} event at time {eventTime}, pos {eventPos}, edge index {edgeIndex}";
				else return $"{eventType} event at time {eventTime}, pos {eventPos}, vertex index {vertexIndex}";
			}
		}

		private PolygonWithHoles inititialPolygon;
		private WavefrontGraph wavefront;
		private StraightSkeleton straightSkeleton;

		private FastStructPriorityQueue<int> eventQueue;

		private int[] wavefrontToStraightSkeletonVertexIndices;

		private float time = 0.0f;
		private readonly float maxEventTime;

		private bool debug;
		private IStraightSkeletonLogger logger = null;


		public StraightSkeletonGenerator(PolygonWithHoles polygonWithHoles, float maxEventTime, bool performDebugChecks)
		{
			this.inititialPolygon = polygonWithHoles;
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

		public void SetLogger(IStraightSkeletonLogger logger)
		{
			this.logger = logger;

			// immediately log the initial input
			logger?.Log($"Initial polygon:", 0);
			logger?.Log($"Outer contour: {string.Join(" ; ", inititialPolygon.outerContourCounterClockwise)}", 1);
			for (int i = 0; i < inititialPolygon.innerContoursClockwise.Count; ++i)
			{
				logger?.Log($"Inner contour #{i}: {string.Join(" ; ", inititialPolygon.innerContoursClockwise[i])}", 1);
			}

			wavefront.SetLogger(logger);
			LogUpcomingEvents();
		}

		public void Step()
		{
			if (!IsDone())
			{
				ProcessEvent();
			}
		}

		public bool IsDone()
		{
			return eventQueue.Count == 0;
		}

		public StraightSkeleton Generate()
		{
			// now we keep popping those events until we're done
			while (eventQueue.Count > 0)
			{
				ProcessEvent();
			}

			return straightSkeleton;
		}

		private QueueEvent DequeueNextEvent()
		{
			int index = eventQueue.Dequeue();
			return GetEventForIndex(index, true);
		}

		private QueueEvent GetEventForIndex(int index, bool resetQueueId)
		{
			if (index >= 0)
			{
				int edgeIndex = index;
				ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
				if (!edgeEvent.eventType.IsValid()) return QueueEvent.Invalid;
				if (resetQueueId) edgeEvent.queueId = -1;
				return new QueueEvent(edgeIndex, edgeEvent);
			}
			else
			{
				int vertexIndex = QueueEvent.QueueIndexToVertexIndex(index);
				ref VertexData vertexData = ref wavefront.vertexDatas[vertexIndex];
				if (!vertexData.partOfSplitEvent) return QueueEvent.Invalid;
				if (resetQueueId) vertexData.queueId = -1;
				return new QueueEvent(vertexIndex, vertexData);
			}
		}

		private void ProcessEvent()
		{
			if (wavefront.nVertices > 10000) throw new ArgumentException($"I think we stumbled upon an infinite loop.");

			// find the first event in the queue
			QueueEvent queueEvent = DequeueNextEvent();
			while (queueEvent.eventType == EventType.None && eventQueue.Count > 0) queueEvent = DequeueNextEvent();
			if (queueEvent.eventType == EventType.None) return;

			// if this is the first event that is over time, we just go over all remaining edges in the queue and spawn them at their position at maxEventTime
			if (queueEvent.eventTime > maxEventTime + Geometry2D.EPS)
			{
				SpawnRemainingEdgesAtMaxTime(queueEvent);
				return;
			}

			logger?.Log($"Process event {queueEvent}", 0);
			switch (queueEvent.eventType)
			{
				case EventType.Edge:
					ProcessEdgeEvent(queueEvent);
					break;

				case EventType.VertexSplit:
					ProcessSplitEvent(queueEvent);
					break;

				case EventType.Nook:
					ProcessNookEvent(queueEvent);
					break;
			}

			// if we emptied the entire set of events that happened at the same time, we update the wavefront
			time = queueEvent.eventTime;
			wavefront.UpdateWavefront(time);

			// validate the graph to catch bugs early
			wavefront.ValidateState(time);

			LogUpcomingEvents();
		}

		private void LogUpcomingEvents()
		{
			// we now log the upcoming events
			if (logger != null)
			{
				int n = math.min(5, eventQueue.Count);
				n = eventQueue.Count;
				int idx = 1;
				logger?.Log($"Upcoming {n} (out of {eventQueue.Count})  events (max allowed {eventQueue.MaxSize}):", 0);
				foreach (int upcomingQueueEventIndex in eventQueue.EnumerateFirst(n))
				{
					var upcomingEvent = GetEventForIndex(upcomingQueueEventIndex, false);
					logger?.Log($"#{idx}: {upcomingEvent}", 1);
					++idx;
				}
			}
		}

		public void AddOrUpdateEdgeEvent(int edgeIndex)
		{

			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			AddOrUpdateEvent(new QueueEvent(edgeIndex, edgeEvent), edgeEvent.queueId, edgeEvent.eventType.IsValid(), out int newQueueId);
			edgeEvent.queueId = newQueueId;
		}

		public void AddOrUpdateVertexEvent(int vertexIndex)
		{

			ref VertexData vertexData = ref wavefront.vertexDatas[vertexIndex];
			bool valid = (vertexData.inWavefront && vertexData.partOfSplitEvent);
			AddOrUpdateEvent(new QueueEvent(vertexIndex, vertexData), vertexData.queueId, valid, out int newQueueId);
			vertexData.queueId = newQueueId;
		}

		private void AddOrUpdateEvent(QueueEvent newEvent, int queueId, bool validEvent, out int newQueueId)
		{
			Profiler.BeginSample("EventQueue.AddOrUpdateEdgeEvent");

			// could be removed and in the batch to be processed at this point
			if (eventQueue.Contains(queueId))
			{
				// if this is a none-event, we either remove it or ignore it
				if (!validEvent)
				{
					logger?.Log($"Existing event in queue {newEvent} got invalidated, so we remove it from the queue", 2);
					eventQueue.Remove(queueId);
					newQueueId = -1;
				}

				// update the priority of the valid event
				else
				{
					logger?.Log($"Existing event in queue {newEvent} got updated to a new event time, so we update the queue", 2);
					eventQueue.UpdatePriority(queueId, newEvent.eventTime);
					newQueueId = queueId;
				}
			}
			else
			{
				if (validEvent)
				{
					logger?.Log($"New event {newEvent} will be added to queue", 2);
					eventQueue.Enqueue(newEvent.index, newEvent.eventTime, out newQueueId);
				}
				else
				{
					//logger?.Log($"New event {newEvent} is invalid, so we don't add it", 0);
					newQueueId = -1;
				}
			}

			Profiler.EndSample();
		}

		private void ProcessEdgeEvent(in QueueEvent queueEvent)
		{
			// we are going to be merging all vertices generated in the wavefront, because their difference don't matter to us
			int edgeIndex = queueEvent.edgeIndex;

			// first add the new vertex to the wavefront and update the wavefront vertices
			int newVertexIndex = wavefront.AddVertexToWavefrontAndRemoveEdge(edgeIndex, true);

			// invalid state - we didn't remove an edge
			if (newVertexIndex == -1) return;

			float2 newVertex = wavefront.vertices[newVertexIndex];
			int newSKVertexIndex = straightSkeleton.AddVertex(newVertex, wavefront.GetVertexTime(newVertexIndex));

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
			if (prevSKVertex.x != newVertex.x || prevSKVertex.y != newVertex.y)
			{

				AddStraightSkeletonArc(prevSKVertexIndex, newSKVertexIndex);
			}
			if (nextSKVertex.x != newVertex.x || nextSKVertex.y != newVertex.y)
			{
				AddStraightSkeletonArc(newSKVertexIndex, nextSKVertexIndex);
			}
		}

		private List<Edge> newStraightSkeletonEdges = new List<Edge>();
		private List<int> splitVertexIndices = new List<int>();
		private List<int> newVertexIndices = new List<int>();
		private void ProcessSplitEvent(in QueueEvent queueEvent)
		{
			// we are going to be merging all vertices generated in the wavefront, because their difference don't matter to us
			newStraightSkeletonEdges.Clear();
			splitVertexIndices.Clear();
			newVertexIndices.Clear();

			// first, make sure this is still a valid spit - another split or event might have deleted some of the edges or involved vertices
			int reflexVertexIndex = queueEvent.vertexIndex;
			if (!wavefront.vertexDatas[reflexVertexIndex].partOfSplitEvent)
			{
				return;
			}

			ref var reflexVertexData = ref wavefront.vertexDatas[reflexVertexIndex];
			int edgeIndex = reflexVertexData.splitEdge;

			float splitTime = reflexVertexData.splitTime;
			float2 splitPos = reflexVertexData.splitPos;

			// add the reflex vertex to the list
			if (reflexVertexData.type != WavefrontVertexType.Reflex) throw new ArgumentException($"Split edge event {reflexVertexData} was triggered by a non-reflex vertex.");
			splitVertexIndices.Add(reflexVertexIndex);

			// now if we're splitting an edge, we first split the edge and add the split point to the list
			// ... but never split the same edge twice!
			if (reflexVertexData.splitPoint == SplitPoint.Edge)
			{
				if (wavefront.edgeEvents[edgeIndex].eventType != EventType.NotInWavefront)
				{
					int newVertexIndex = wavefront.SplitEdge(reflexVertexIndex, edgeIndex, splitTime, splitPos);
					splitVertexIndices.Add(newVertexIndex);
					logger?.Log($"Split point is {reflexVertexData.splitPoint}, so we first split edge {edgeIndex} and add a new vertex {newVertexIndex} at the split point, which is now part of the set of split vertex indices, to be removed later.", 2);
				}
			}

			// if we're not splitting an edge for we have a split at one of the edge ends, we actually add the endpoint as well if it is NOT a relfex vertex
			// this is because when two reflex vertices collide, they will also register to each other as a split point, but will already have been added earlier
			else
			{
				ref var edge = ref wavefront.edges[edgeIndex];
				int splitVertexIndex = edge.GetPoint(reflexVertexData.splitPoint);
				ref var splitVertexData = ref wavefront.vertexDatas[splitVertexIndex];
				if (splitVertexData.type == WavefrontVertexType.Convex && !splitVertexIndices.Contains(splitVertexIndex))
				{
					logger?.Log($"Split point is {reflexVertexData.splitPoint} of edge {edgeIndex}, so we intersect exactly at existing vertex {splitVertexIndex}, which is now part of the set of split vertex indices, to be removed later.", 2);
					splitVertexIndices.Add(splitVertexIndex);
				}
			}

			// go over all vertices and see if they split at the same point
			for (int i = 0; i < wavefront.nVertices; ++i)
			{
				if (i == reflexVertexIndex) continue;

				ref var otherVertexData = ref wavefront.vertexDatas[i];

				if (otherVertexData.partOfSplitEvent && Geometry2D.CompareToLowPrecision(reflexVertexData.splitPos, otherVertexData.splitPos) == 0)
				{
					splitVertexIndices.Add(i);
					if (otherVertexData.queueId != -1)
					{
						eventQueue.Remove(otherVertexData.queueId);
						otherVertexData.queueId = -1;
					}
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
				if (prevSKVertex.x != nextSKVertex.x || prevSKVertex.y != nextSKVertex.y)
				{
					AddStraightSkeletonArc(prevSKVertexIndex, nextSKVertexIndex);

				}
			}
		}

		private void ProcessNookEvent(in QueueEvent queueEvent)
		{
			// get the event and do a sanity check
			int edgeIndex = queueEvent.edgeIndex;
			ref var edgeEvent = ref wavefront.edgeEvents[edgeIndex];
			ref var edge = ref wavefront.edges[edgeIndex];
			if (edgeEvent.eventType != EventType.Nook) throw new InvalidOperationException($"There is no nook event associated with edge {edge}.");

			// remove the nook from the wavefront
			wavefront.RemoveNook(edgeIndex);

			// this is easy - just spawn the edge, since the vertices already exist
			int prevSKVertexIndex = GetStraightSkeletonVertexAtTime(edge.prevVertexIndex, edgeEvent.eventTime);
			int nextSKVertexIndex = GetStraightSkeletonVertexAtTime(edge.nextVertexIndex, edgeEvent.eventTime);

			// see if this warrants an edge
			float2 prevSKVertex = straightSkeleton.vertices[prevSKVertexIndex];
			float2 nextSKVertex = straightSkeleton.vertices[nextSKVertexIndex];

			// add the edges to the straight skeleton
			if (prevSKVertex.x != nextSKVertex.x || prevSKVertex.y != nextSKVertex.y)
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

			logger?.Log($"Reached maximum event time {maxEventTime}, spawning remaining edges at their position at that time...", 0);

			// we generate straight skeleton edges from the current vertex to its position at the max time,
			// and also update the mapping from the old SK vertex position to the new one
			for (int i = 0; i < wavefront.nVertices; ++i)
			{
				ref VertexData vertexData = ref wavefront.vertexDatas[i];
				if (!vertexData.inWavefront) continue;
				SpawnEdgeFromPosToMaxTimePos(i);
			}

			// re-use this list for convenience
			for (int i = 0; i < wavefront.edgeEvents.Length; ++i)
			{
				ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[i];
				if (edgeEvent.eventType != EventType.NotInWavefront)
				{
					SpawnEdgeAtMaxTime(i);
				}
			}

			// clear the event queue - we're done!
			eventQueue.Clear();
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
			if (vertexData.creationTime < maxEventTime - Geometry2D.EPS)
			{
				int startSKVertexIndex = GetStraightSkeletonVertexAtCreationTime(vertexIndex);
				wavefrontToStraightSkeletonVertexIndices[vertexIndex] = -1;
				int endSKVertexIndex = GetStraightSkeletonVertexAtMaxTime(vertexIndex);
				AddStraightSkeletonArc(startSKVertexIndex, endSKVertexIndex);
			}
		}

		public void DrawGizmos(Camera camera)
		{
			Gizmos.color = Color.red;

			for (int i = 0; i < straightSkeleton.nEdges; ++i)
			{
				var edge = straightSkeleton.edges[i];
				var prevVertex = straightSkeleton.vertices[edge.prevVertexIndex];
				var nextVertex = straightSkeleton.vertices[edge.nextVertexIndex];

				Gizmos.DrawLine(new Vector3(prevVertex.x, straightSkeleton.vertexTimes[edge.prevVertexIndex], prevVertex.y), new Vector3(nextVertex.x, straightSkeleton.vertexTimes[edge.nextVertexIndex], nextVertex.y));
			}

			if (!IsDone()) wavefront.DrawGizmos(time, camera);
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