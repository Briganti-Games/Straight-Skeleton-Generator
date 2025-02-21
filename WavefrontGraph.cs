using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Profiling;
using UnityEngine.Rendering;
using UnityEngine.UI;

namespace Briganti.StraightSkeletonGeneration
{
	public interface IEventTimeChangeListener
	{
		void AddOrUpdateEdgeEvent(int edgeIndex);
	}

	public class WavefrontGraph : VertexGraph
	{
		public readonly EdgeEvent[] edgeEvents;
		public readonly VertexData[] vertexDatas;
		private readonly IEventTimeChangeListener eventTimeListener;

		private HashSet<int> affectedEdges = new HashSet<int>();
		private HashSet<int> affectedVertices = new HashSet<int>();

		private bool debug;


		public WavefrontGraph(PolygonWithHoles polygonWithHoles, IEventTimeChangeListener eventTimeListener, bool debug)
		{
			this.debug = debug;
			this.eventTimeListener = eventTimeListener;

			int nVertices = polygonWithHoles.outerContour.Length;
			for (int i = 0; i < polygonWithHoles.innerContours.Count; ++i)
			{
				nVertices += polygonWithHoles.innerContours[i].Length;
			}

			int nEdges = nVertices;

			// we now add all the edges that will be added for the straight skeleton - the exact number is known
			// each vertex can spawn at most another vertex, because it either collapses an edge, or is moved when it only collapses after the max event time
			int maxVertices = nVertices * 2;
			int nArcs = ((2 * nVertices) - 3);
			int maxEdges = nEdges + nArcs;

			// this can be improved and the upper limit can be calculated EXACTLY - see the lemmas in Stefan Huber's PhD
			Initialize(maxVertices, maxEdges);

			// we store additional data for each edge and vertex
			edgeEvents = new EdgeEvent[maxEdges];
			vertexDatas = new VertexData[maxVertices];

			// add all contours
			AddContour(polygonWithHoles.outerContour);
			for (int i = 0; i < polygonWithHoles.innerContours.Count; ++i)
			{
				AddContour(polygonWithHoles.innerContours[i]);
			}
		}

		public void GenerateInitialEvents()
		{

			// update all events
			for (int i = 0; i < nEdges; ++i)
			{
				UpdateEventTime(i, 0);
			}
		}

		private void AddContour(float2[] contour)
		{
			int vertexOffset = nVertices;
			int edgeOffset = nEdges;
			for (int i = 0; i < contour.Length; ++i)
			{
				float2 curr = contour[i];

				// add the vertex
				AddVertex(curr);

				int prevVertexIndex = vertexOffset + (i - 1 + contour.Length) % contour.Length;
				int nextVertexIndex = vertexOffset + (i + 1) % contour.Length;

				// set the vertex data
				int prevEdgeIndex = edgeOffset + (i - 1 + contour.Length) % contour.Length;
				int nextEdgeIndex = edgeOffset + i;
				vertexDatas[i].nextSplitReflexVertexIndex = -1; // we chain all reflex vertices that split an edge at the same position
				vertexDatas[i].UpdateConnections(prevVertexIndex, nextVertexIndex, prevEdgeIndex, nextEdgeIndex);

				// also store the edge that connects this vertex and the next one
				int edgeIndex = AddEdge(new Edge(i, (i + 1) % contour.Length));

				// by default there is no event for this edge
				edgeEvents[edgeIndex].eventType = EventType.None;
			}

			// now update the vertex data - calculate velocity etc now that the entire edge is known
			for (int i = 0; i < contour.Length; ++i)
			{
				UpdateVertexData(vertexOffset + i, 0);
			}
		}

		public int AddVertexToWavefrontAndRemoveEdge(int disappearedEdgeIndex, bool updateEdgeEvents)
		{
			// get the edge event
			ref EdgeEvent edgeEvent = ref edgeEvents[disappearedEdgeIndex];

			// if this event has changed by now, we ignore it
			if (edgeEvent.eventType != EventType.Edge) return -1;

			float2 pos = edgeEvent.eventPos;
			float time = edgeEvent.eventTime;

			ref Edge disappearedEdge = ref edges[disappearedEdgeIndex];
			ref VertexData prevVertexData = ref vertexDatas[disappearedEdge.prevVertexIndex];
			ref VertexData nextVertexData = ref vertexDatas[disappearedEdge.nextVertexIndex];

			// this is a special case - we are at the end of the wavefront and it's been reduced to exactly one edge (a nook)
			// in this case, we skip this event entirely, it will be handled later
			if (prevVertexData.prevVertexIndex == disappearedEdge.nextVertexIndex && nextVertexData.nextVertexIndex == disappearedEdge.prevVertexIndex)
			{
				return -1;
			}


			if (prevVertexData.nextEdgeIndex != disappearedEdgeIndex) throw new ArgumentException($"Vertex {disappearedEdge.prevVertexIndex} is not properly linked to edge {disappearedEdgeIndex}");
			if (nextVertexData.prevEdgeIndex != disappearedEdgeIndex) throw new ArgumentException($"Vertex {disappearedEdge.nextVertexIndex} is not properly linked to edge {disappearedEdgeIndex}");

			// we flag that these vertices were removed from the wavefront, since they merged into a new vertex
			prevVertexData.inWavefront = false;
			nextVertexData.inWavefront = false;

			int prevEdgeIndex = prevVertexData.prevEdgeIndex;
			int nextEdgeIndex = nextVertexData.nextEdgeIndex;

			int prevVertexIndex = prevVertexData.prevVertexIndex;
			int nextVertexIndex = nextVertexData.nextVertexIndex;

			ref Edge prevEdge = ref edges[prevVertexData.prevEdgeIndex];
			ref Edge nextEdge = ref edges[nextVertexData.nextEdgeIndex];

			// firstly, just add the vertex
			int newVertexIndex = AddVertex(pos);

			// update the new vertex with the right connections
			UpdateConnections(newVertexIndex, prevVertexIndex, nextVertexIndex, prevEdgeIndex, nextEdgeIndex);
			vertexDatas[newVertexIndex].creationTime = time;

			return newVertexIndex;
		}

		public void RemoveNook(int edgeIndex)
		{
			ref EdgeEvent edgeEvent = ref edgeEvents[edgeIndex];
			if (edgeEvent.eventType != EventType.Nook) throw new ArgumentException($"Edge event {edgeEvent} is not a nook event!");

			ref Edge edge = ref edges[edgeIndex];

			// remove the edge from the wavefront (and along it also the other edge that forms the same nook)
			vertexDatas[edge.prevVertexIndex].inWavefront = false;
			vertexDatas[edge.nextVertexIndex].inWavefront = false;
		}

		public void SplitEdge(int splitEdgeIndex, List<int> newVertexIndices)
		{
			// get the edge event
			ref EdgeEvent edgeEvent = ref edgeEvents[splitEdgeIndex];

			// if this event has changed by now, we ignore it
			if (edgeEvent.eventType != EventType.Split) return;

			float2 pos = edgeEvent.eventPos;
			float time = edgeEvent.eventTime;

			ref Edge oldEdge = ref edges[splitEdgeIndex];

			// we let the reflex vertex die - it served its purpose
			int reflexVertexIndex = edgeEvent.firstSplitReflexVertexIndex;
			ref VertexData reflexVertexData = ref vertexDatas[reflexVertexIndex];

			// if the reflex vertex was already removed from the wavefront, this is the second split event for the same reflex vertex at the same time
			// this can only mean that the reflex vertex hit 2 edges at their meeting point, and we already processed the split
			if (!reflexVertexData.inWavefront) return;

			reflexVertexData.inWavefront = false;
			if (reflexVertexData.nextSplitReflexVertexIndex != -1) throw new NotSupportedException($"Multi-split events are currently not yet supported.");

			// create two new vertices - one to be part of the "prev" polygon, one to be part of the "next" polygon
			int newPrevVertexIndex = AddVertex(pos);
			int newNextVertexIndex = AddVertex(pos);
			vertexDatas[newPrevVertexIndex].creationTime = time;
			vertexDatas[newNextVertexIndex].creationTime = time;

			// might be more later for multisplit events
			newVertexIndices.Add(newPrevVertexIndex);
			newVertexIndices.Add(newNextVertexIndex);

			// we create two new edges from the old edge
			int newPrevEdgeIndex = AddEdge(oldEdge.prevVertexIndex, newPrevVertexIndex);
			int newNextEdgeIndex = AddEdge(newNextVertexIndex, oldEdge.nextVertexIndex);

			// we now need to update the vertex connections to and from the split vertex point
			UpdateConnections(newPrevVertexIndex, oldEdge.prevVertexIndex, reflexVertexData.nextVertexIndex, newPrevEdgeIndex, reflexVertexData.nextEdgeIndex);
			UpdateConnections(newNextVertexIndex, reflexVertexData.prevVertexIndex, oldEdge.nextVertexIndex, reflexVertexData.prevEdgeIndex, newNextEdgeIndex);
		}

		public int SplitGraphAtVertex(int splitEdgeIndex, List<int> newVertexIndices)
		{
			// get the edge event
			ref EdgeEvent edgeEvent = ref edgeEvents[splitEdgeIndex];

			// if this event has changed by now, we ignore it
			if (edgeEvent.eventType != EventType.SplitAtVertex) return -1;

			float2 pos = edgeEvent.eventPos;
			float time = edgeEvent.eventTime;

			ref Edge oldEdge = ref edges[splitEdgeIndex];

			// we let the reflex vertex die - it served its purpose
			int reflexVertexIndex = edgeEvent.firstSplitReflexVertexIndex;
			ref VertexData reflexVertexData = ref vertexDatas[reflexVertexIndex];

			// if the reflex vertex was already removed from the wavefront, this is the second split event for the same reflex vertex at the same time
			// this can only mean that the reflex vertex hit 2 edges at their meeting point, and we already processed the split
			if (!reflexVertexData.inWavefront) return -1;

			reflexVertexData.inWavefront = false;
			if (reflexVertexData.nextSplitReflexVertexIndex != -1) throw new NotSupportedException($"Multi-split events are currently not yet supported.");

			// see which vertex we're splitting the graph over
			float2 oldEdgePrevVertex = GetVertexPosAtTime(oldEdge.prevVertexIndex, time);
			float2 oldEdgeNextVertex = GetVertexPosAtTime(oldEdge.nextVertexIndex, time);

			int splitVertexIndex = -1;
			if (math.distancesq(oldEdgePrevVertex, pos) < Geometry.EPSSQ) splitVertexIndex = oldEdge.prevVertexIndex;
			else if (math.distancesq(oldEdgeNextVertex, pos) < Geometry.EPSSQ) splitVertexIndex = oldEdge.nextVertexIndex;
			else throw new ArgumentException($"This is a vertex split event, but the reflex position {pos} is at neither endpoints {oldEdgeNextVertex} or {oldEdgeNextVertex} of edge {oldEdge}");
			ref VertexData splitVertexData = ref vertexDatas[splitVertexIndex];
			splitVertexData.inWavefront = false;

			// create two new vertices - one to be part of the "prev" polygon, one to be part of the "next" polygon
			int newPrevVertexIndex = AddVertex(pos);
			int newNextVertexIndex = AddVertex(pos);
			vertexDatas[newPrevVertexIndex].creationTime = time;
			vertexDatas[newNextVertexIndex].creationTime = time;

			// might be more later for multisplit events
			newVertexIndices.Add(newPrevVertexIndex);
			newVertexIndices.Add(newNextVertexIndex);

			// we now need to update the vertex connections to and from the split vertex point
			UpdateConnections(newPrevVertexIndex, splitVertexData.prevVertexIndex, reflexVertexData.nextVertexIndex, splitVertexData.prevEdgeIndex, reflexVertexData.nextEdgeIndex);
			UpdateConnections(newNextVertexIndex, reflexVertexData.prevVertexIndex, splitVertexData.nextVertexIndex, reflexVertexData.prevEdgeIndex, splitVertexData.nextEdgeIndex);

			return splitVertexIndex;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private void UpdateConnections(int vertexIndex, int prevVertexIndex, int nextVertexIndex, int prevEdgeIndex, int nextEdgeIndex)
		{
			// update the vertex itself
			ref var vertexData = ref vertexDatas[vertexIndex];
			vertexData.UpdateConnections(prevVertexIndex, nextVertexIndex, prevEdgeIndex, nextEdgeIndex);

			// also update the adjacent edges
			edges[vertexData.prevEdgeIndex].nextVertexIndex = vertexIndex;
			edges[vertexData.nextEdgeIndex].prevVertexIndex = vertexIndex;

			// also update the adjacent vertices
			vertexDatas[prevVertexIndex].nextVertexIndex = vertexIndex;
			vertexDatas[nextVertexIndex].prevVertexIndex = vertexIndex;
			vertexDatas[prevVertexIndex].nextEdgeIndex = prevEdgeIndex;
			vertexDatas[nextVertexIndex].prevEdgeIndex = nextEdgeIndex;

			// remember these are affected by graph changes, so we can update them later
			affectedVertices.Add(vertexIndex);
			affectedVertices.Add(prevVertexIndex);
			affectedVertices.Add(nextVertexIndex);

			affectedEdges.Add(prevEdgeIndex);
			affectedEdges.Add(nextEdgeIndex);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private bool IsEdgeInWavefront(int edgeIndex) {
			ref var edge = ref edges[edgeIndex];
			return vertexDatas[edge.prevVertexIndex].inWavefront && vertexDatas[edge.nextVertexIndex].inWavefront;
		}

		public void UpdateWavefront(float time)
		{
			foreach (int vertexIndex in affectedVertices)
			{
				if (vertexDatas[vertexIndex].inWavefront)
				{
					UpdateVertexData(vertexIndex, time);
				}
			}
			affectedVertices.Clear();

			foreach (int edgeIndex in affectedEdges)
			{
				UpdateEventTime(edgeIndex, time);
			}
			affectedEdges.Clear();
		}

		public int GetNextVertexIndex(int vertexIndex)
		{
			ref VertexData vertexData = ref vertexDatas[vertexIndex];
			ref Edge edge = ref edges[vertexData.nextEdgeIndex];
			return edge.nextVertexIndex;
		}

		public int GetPrevVertexIndex(int vertexIndex)
		{
			ref VertexData vertexData = ref vertexDatas[vertexIndex];
			ref Edge edge = ref edges[vertexData.prevEdgeIndex];
			return edge.prevVertexIndex;
		}

		private void UpdateVertexData(int vertexIndex, float time)
		{

			float2 vertex = GetVertexPosAtTime(vertexIndex, time);
			ref VertexData vertexData = ref vertexDatas[vertexIndex];
			float2 prevVertex = GetVertexPosAtTime(vertexData.prevVertexIndex, time);
			float2 nextVertex = GetVertexPosAtTime(vertexData.nextVertexIndex, time);

			float2 velocity = CalculateVelocity(vertexIndex, prevVertex, vertex, nextVertex);
			bool isReflex = Geometry.IsRelfexVertex(prevVertex, vertex, nextVertex);
			var type = (isReflex ? WavefrontVertexType.Reflex : WavefrontVertexType.Convex);

			// the vertex lies on a parallel line - we can't calculate a real velocity
			if (Geometry.IsParallelLines(prevVertex, vertex, prevVertex, nextVertex))
			{
				velocity = float2.zero;
			}

			vertexData.velocity = velocity;
			vertexData.type = type;
		}

		public void UpdateEventTime(int edgeIndex, float time)
		{
			// default to max
			ref Edge edge = ref edges[edgeIndex];
			ref EdgeEvent eventData = ref edgeEvents[edgeIndex];

			eventData.Reset();

			// the edge must be either 100% in the current wavefront or removed from it
			bool vertexInWavefront1 = vertexDatas[edge.prevVertexIndex].inWavefront;
			bool vertexInWavefront2 = vertexDatas[edge.nextVertexIndex].inWavefront;
			if (vertexInWavefront1 != vertexInWavefront2) throw new InvalidOperationException($"Edge {edge} is partly in the wavefront and partly not, because it connects vertices {vertexDatas[edge.prevVertexIndex]} and {vertexDatas[edge.nextVertexIndex]}. This is illegal!");
			if (!vertexInWavefront1) return;

			// first, see if we're at the end if a wavefront and it has collapsed to one edge - a nook event
			bool wasNook = UpdateNookEventTime(ref edge, ref eventData, time);
			if (!wasNook)
			{

				// then, the easy part - calculate the edge event
				UpdateEdgeEventTime(ref edge, ref eventData, time);

				// then we see if this edge is being split by a vertex
				UpdateSplitEventTime(edgeIndex, ref eventData, time);
			}

			// add to the queue
			eventTimeListener.AddOrUpdateEdgeEvent(edgeIndex);
		}

		private bool UpdateNookEventTime(ref Edge edge, ref EdgeEvent eventData, float time)
		{
			ref VertexData prevVertexData = ref vertexDatas[edge.prevVertexIndex];
			ref VertexData nextVertexData = ref vertexDatas[edge.nextVertexIndex];

			// this is a special case - we are at the end of the wavefront and it's been reduced to exactly one edge (a nook)
			// in this case, we skip this event entirely, it will be handled later
			if (prevVertexData.prevVertexIndex == edge.nextVertexIndex && nextVertexData.nextVertexIndex == edge.prevVertexIndex)
			{
				// because a nook event always comes in pairs, we check whether the OTHER edge has already been flagged as a nook
				// in that wase, we just skip this one
				ref var otherEdgeEvent = ref edgeEvents[prevVertexData.prevEdgeIndex];
				if (otherEdgeEvent.eventType == EventType.Nook)
				{
					eventData.eventType = EventType.None;
					eventData.eventTime = float.MaxValue;
					eventData.eventPos = float2.zero;
				}

				// first of the two edges we encountered - it's a nook!
				else
				{
					eventData.eventType = EventType.Nook;
					eventData.eventTime = time;
					eventData.eventPos = float2.zero;
				}

				return true;
			}

			return false;
		}

		private void UpdateEdgeEventTime(ref Edge edge, ref EdgeEvent eventData, float _)
		{
			// we store the edge event time in the prevVertex associated with the edge!
			float2 prevVertex = vertices[edge.prevVertexIndex];
			float2 nextVertex = vertices[edge.nextVertexIndex];
			ref VertexData prevData = ref vertexDatas[edge.prevVertexIndex];
			ref VertexData nextData = ref vertexDatas[edge.nextVertexIndex];

			if (Geometry.GetLineIntersection(prevVertex, prevVertex + prevData.velocity, nextVertex, nextVertex + nextData.velocity, out float t0, out float t1))
			{
				if (t0 > -Geometry.EPS && t1 > -Geometry.EPS)
				{
					eventData.eventType = EventType.Edge;
					eventData.eventPos = prevVertex + prevData.velocity * t0;

					float creationTime = Mathf.Max(prevData.creationTime, nextData.creationTime);

					// we "fast forward" the edge into the current timeframe, so we can see where it is and properly calculate
					// the REMAINING time it takes to geft to the collapse point.
					float2 prevVertexCurrentPos = GetVertexPosAtTime(prevVertex, prevData, creationTime);
					float2 nextVertexCurrentPos = GetVertexPosAtTime(nextVertex, nextData, creationTime);

					float2 projPoint = Geometry.ProjectPointOnLine(prevVertexCurrentPos, nextVertexCurrentPos, eventData.eventPos, out float t);

					// the event time is the time it takes for the wavefront to reach this position, AFTER both vertices of the edge were spawned
					eventData.eventTime = creationTime + math.distance(projPoint, eventData.eventPos);
				}
			}
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public float2 GetVertexPosAtTime(int vertexIndex, float time)
		{
			float2 vertex = vertices[vertexIndex];
			ref VertexData vertexData = ref vertexDatas[vertexIndex];
			return GetVertexPosAtTime(vertex, vertexData, time);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private float2 GetVertexPosAtTime(in float2 vertex, in VertexData vertexData, float time)
		{
			//if (time < vertexData.creationTime - Geometry.EPS) throw new ArgumentException($"Time {time} lies in the past of the creation point of vertex {vertex} ({vertexData.creationTime}), this should not be possible.");
			if (time < vertexData.creationTime - Geometry.EPS) time = vertexData.creationTime;
			time -= vertexData.creationTime;
			return vertex + vertexData.velocity * time;
		}

		private void UpdateSplitEventTime(int edgeIndex, ref EdgeEvent edgeEvent, float time)
		{
			// go over all active reflex vertices and see if they split us
			for (int i = 0; i < nVertices; ++i)
			{
				if (vertexDatas[i].inWavefront && vertexDatas[i].type == WavefrontVertexType.Reflex)
				{
					UpdateSplitEventTime(edgeIndex, i, ref edgeEvent, time);
				}
			}
		}

		private void UpdateSplitEventTime(int edgeIndex, int reflexVertexIndex, ref EdgeEvent edgeEvent, float currentTime)
		{
			ref Edge edge = ref edges[edgeIndex];

			// can't be adjacent
			if (edge.prevVertexIndex == reflexVertexIndex || edge.nextVertexIndex == reflexVertexIndex) return;

			float2 prevVertex = GetVertexPosAtTime(edge.prevVertexIndex, currentTime);
			float2 nextVertex = GetVertexPosAtTime(edge.nextVertexIndex, currentTime);

			float2 reflexVertex = GetVertexPosAtTime(reflexVertexIndex, currentTime);
			ref VertexData reflexVertexData = ref vertexDatas[reflexVertexIndex];
			float2 reflexVelocity = reflexVertexData.velocity;

			// if the reflex vertex lies ON the prev or next vertex, we don't have a split
			if (math.distancesq(prevVertex, reflexVertex) < Geometry.EPSSQ || math.distancesq(nextVertex, reflexVertex) < Geometry.EPSSQ) return;

			// calculate the normal of the edge
			float2 normal = Geometry.Rotate90DegreesClockwise(nextVertex - prevVertex);

			// calculate the angle between the normal and reflex vertex
			float angle = Geometry.Angle(-normal, reflexVelocity);

			// if the reflex vertex is moving in the opposite direction, the reflex vertex will be swallowed by the wavefront before a split can occur with this edge
			if (angle > Mathf.PI / 2) return;

			// we can now "easily" calculate the distance the edge needs to travel before it hits our reflex vertex
			// (it took 6 hours to find the solution to this problem)
			float2 reflexVertexProj = Geometry.ProjectPointOnLine(prevVertex, nextVertex, reflexVertex, out _);
			float H = math.distance(reflexVertex, reflexVertexProj);
			float V = math.length(reflexVelocity);

			// this is now the distance the edge travelled along its normal to reach the collision point
			float L = H / (1 + V * math.cos(angle));
			float time = L;

			// we can now determine the position where the edge meets the reflex vertex
			//float distanceAlongVertexLine = L * V;
			float2 eventPos = reflexVertex + reflexVelocity * time;

			// because lines always move at constant speed, the time is identical to the move distance of the edge
			float eventTime = currentTime + time;

			// there was another event earlier
			if (edgeEvent.eventType != EventType.None && edgeEvent.eventTime < eventTime) return;

			// if the reflex vertex collides with us after our collapse, we don't bother
			if (edgeEvent.eventType == EventType.Edge && eventTime > edgeEvent.eventTime - Geometry.EPS) return;

			// see if we are in the same time & place for the split
			bool sameTimeAndPlace = (math.abs(eventTime - edgeEvent.eventTime) < Geometry.EPS) && (math.distancesq(eventPos, edgeEvent.eventPos) < Geometry.EPSSQ);

			// we have multiple splits with this edge at the same time - we just skip this one and do them one by one
			if (edgeEvent.eventType == EventType.Split && !sameTimeAndPlace) return;

			// now finally, we need to make sure that this still falls into the edge once it moved that far
			float2 prevVertexAtEventTime = GetVertexPosAtTime(edge.prevVertexIndex, eventTime);
			float2 nextVertexAtEventTime = GetVertexPosAtTime(edge.nextVertexIndex, eventTime);
			Geometry.ProjectPointOnLine(prevVertexAtEventTime, nextVertexAtEventTime, eventPos, out float tLine);
			if (tLine < 0 || tLine > 1) return;

			// at this point, if we have a split, it MUST be a multi split!
			if (edgeEvent.eventType == EventType.Split && !sameTimeAndPlace) throw new ArgumentException($"We should have a split event, but current even {edgeEvent} does not occur at time {eventTime} or place {eventPos}.");

			// we have a split event!
			edgeEvent.eventType = EventType.Split;
			edgeEvent.eventTime = eventTime;
			edgeEvent.eventPos = eventPos;

			// if we split at exactly the corner, we have a special type of split event - this event will NOT generate new edges,
			// but will still split the entire graph into two pieces and spawn new vertices at the split point.
			if (tLine < Geometry.EPS || tLine > 1.0f - Geometry.EPS) edgeEvent.eventType = EventType.SplitAtVertex;

			// this is the first split event at this time & place, store the index
			if (edgeEvent.firstSplitReflexVertexIndex == -1)
			{
				edgeEvent.firstSplitReflexVertexIndex = reflexVertexIndex;
			}
			else
			{
				AddSplitReflexIndexToChain(ref edgeEvent, reflexVertexIndex);
			}
		}

		private void AddSplitReflexIndexToChain(ref EdgeEvent edgeEvent, int reflexVertexIndex)
		{
			ref VertexData vertexData = ref vertexDatas[edgeEvent.firstSplitReflexVertexIndex];

			while (vertexData.nextSplitReflexVertexIndex != -1)
			{
				vertexData = ref vertexDatas[vertexData.nextSplitReflexVertexIndex];
			}

			vertexData.nextSplitReflexVertexIndex = reflexVertexIndex;
		}

		public static float2 CalculateVelocity(int vertexIndex, in float2 prev, in float2 curr, in float2 next)
		{
			float2 prevToCurr = curr - prev;
			float2 nextToCurr = next - curr;

			// rotate them 90°
			float2 moveDir1 = math.normalize(Geometry.RotateMinus90Degrees(prevToCurr));
			float2 moveDir2 = math.normalize(Geometry.RotateMinus90Degrees(nextToCurr));

			float2 dir = math.normalize(moveDir1 + moveDir2);
			float speed = 2.0f / (math.dot(moveDir1, dir) + math.dot(moveDir2, dir));

			float2 velocity = dir * speed;
			//Debug.Log("Vertex #" + vertexIndex + " has dir " + dir + " and speed " + speed + " based on " + moveDir1 + " and " + moveDir2 + " and final dir has angle " + (Geometry.GetAngle(dir) * Mathf.Rad2Deg) + ", angle between two adjacent dirs is " + (Geometry.SignedAngle(prevToCurr, nextToCurr) * Mathf.Rad2Deg));
			return velocity;
		}

		public float GetVertexTime(int vertexIndex)
		{
			return vertexDatas[vertexIndex].creationTime;
		}

		public void ValidateState(float time)
		{
			if (!debug) return;
			for (int i = 0; i < nVertices; ++i)
			{
				ref var vertexData = ref vertexDatas[i];
				if (!vertexData.inWavefront) continue;

				ref var prevEdge = ref edges[vertexData.prevEdgeIndex];
				ref var nextEdge = ref edges[vertexData.nextEdgeIndex];

				if (prevEdge.nextVertexIndex != i) throw new ArgumentException($"Prev edge {prevEdge} (idx {vertexData.prevEdgeIndex}) does not connect correctly with vertex {vertexData} (idx {i})");
				if (nextEdge.prevVertexIndex != i) throw new ArgumentException($"Next edge {nextEdge} (idx {vertexData.nextEdgeIndex}) does not connect correctly with vertex {vertexData} (idx {i})");
				if (prevEdge.prevVertexIndex != vertexData.prevVertexIndex) throw new ArgumentException($"Prev edge {prevEdge} (idx {vertexData.prevEdgeIndex}) does not have the same previous vertex {prevEdge.prevVertexIndex} defined as the vertex {vertexData} (idx {i})");
				if (nextEdge.nextVertexIndex != vertexData.nextVertexIndex) throw new ArgumentException($"Nexts edge {nextEdge} (idx {vertexData.nextEdgeIndex}) does not have the same previous vertex {prevEdge.nextVertexIndex} defined as the vertex {vertexData} (idx {i})");
			}

			string GetVertexDescription(int vertexIndex)
			{
				return $"{vertexIndex} ({GetVertexPosAtTime(vertexIndex, time)})";
			}

			// print the vertex chains
			StringBuilder ss = new StringBuilder();
			ss.AppendLine("Chains at time " + time + ":");
			int chainCount = 0;
			HashSet<int> visitedVertices = new HashSet<int>();
			for (int i = 0; i < nVertices; ++i)
			{
				ref var vertexData = ref vertexDatas[i];
				if (!vertexData.inWavefront) continue;
				if (visitedVertices.Contains(i)) continue;

				int initialVertexIndex = i;
				int vertexIndex = initialVertexIndex;
				ss.Append("Chain #" + (++chainCount) + ":");
				ss.Append(GetVertexDescription(initialVertexIndex));
				visitedVertices.Add(initialVertexIndex);
				do
				{
					vertexIndex = vertexData.nextVertexIndex;
					vertexData = ref vertexDatas[vertexIndex];
					ss.Append(" -> ");
					ss.Append(GetVertexDescription(vertexIndex));
					visitedVertices.Add(vertexIndex);
				}
				while (vertexIndex != initialVertexIndex);
				ss.AppendLine();
			}

			Debug.Log(ss.ToString());
		}



		public void DrawGizmos(float time)
		{
			Gizmos.color = Color.cyan;

			for (int i = 0; i < nEdges; ++i)
			{
				ref var edge = ref edges[i];
				if (vertexDatas[edge.prevVertexIndex].inWavefront && vertexDatas[edge.prevVertexIndex].nextEdgeIndex == i && vertexDatas[edge.nextVertexIndex].inWavefront && vertexDatas[edge.nextVertexIndex].prevEdgeIndex == i)
				{
					var prevVertex = GetVertexPosAtTime(edge.prevVertexIndex, time);
					var nextVertex = GetVertexPosAtTime(edge.nextVertexIndex, time);
					float prevCreationTime = vertexDatas[edge.prevVertexIndex].creationTime;
					float nextCreationTime = vertexDatas[edge.nextVertexIndex].creationTime;

					Gizmos.DrawLine(new Vector3(prevVertex.x, prevCreationTime, prevVertex.y), new Vector3(nextVertex.x, nextCreationTime, nextVertex.y));
				}
			}

			Gizmos.color = Color.gray;

			for (int i = 0; i < nVertices; ++i)
			{
				var vertex = GetVertexPosAtTime(i, time);
				ref var vertexData = ref vertexDatas[i];
				if (!vertexData.inWavefront) continue;

				var endPoint = vertex + math.normalize(vertexData.velocity) * 0.3f;
				float creationTime = vertexDatas[i].creationTime;

				Gizmos.DrawLine(new Vector3(vertex.x, creationTime + 0.01f, vertex.y), new Vector3(endPoint.x, creationTime + 0.01f, endPoint.y));
			}

		}
	}
}