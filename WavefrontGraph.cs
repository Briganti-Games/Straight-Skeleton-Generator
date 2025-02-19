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
	public interface IEventTimeChangeListener
	{
		void AddEdgeEvent(int edgeIndex);
		void UpdateEdgeEvent(int edgeIndex);
	}

	public class WavefrontGraph : VertexGraph
	{
		public readonly EdgeEvent[] edgeEvents;
		public readonly VertexData[] vertexDatas;
		private readonly IEventTimeChangeListener eventTimeListener;


		public WavefrontGraph(PolygonWithHoles polygonWithHoles, IEventTimeChangeListener eventTimeListener)
		{
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
				eventTimeListener.AddEdgeEvent(i);
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
			vertexDatas[newVertexIndex].UpdateConnections(prevVertexIndex, nextVertexIndex, prevEdgeIndex, nextEdgeIndex);
			vertexDatas[newVertexIndex].creationTime = time;

			// now we cut out the disappeared edge from the graph and immediately connect the adjacent vertices
			prevEdge.nextVertexIndex = newVertexIndex;
			nextEdge.prevVertexIndex = newVertexIndex;

			// if there are two remaining edges, we are at the end of our process here (the graph is only 2 edges long and they connect to each other),
			// we add an additional nook event that connects the final two vertices
			if (prevEdge.nextVertexIndex == nextEdge.prevVertexIndex && prevEdge.prevVertexIndex == nextEdge.nextVertexIndex)
			{

				ref EdgeEvent prevEdgeEventData = ref edgeEvents[prevVertexData.prevEdgeIndex];
				ref EdgeEvent nextEdgeEventData = ref edgeEvents[nextVertexData.nextEdgeIndex];

				float2 p1 = vertices[prevEdge.nextVertexIndex];
				float2 p2 = vertices[prevEdge.prevVertexIndex];

				// disable these events
				prevEdgeEventData.Reset();
				nextEdgeEventData.Reset();

				// ... but we might have to create one final nook edge though!
				if (math.distance(p1, p2) > Geometry.EPS)
				{
					prevEdgeEventData.eventType = EventType.Nook;
					prevEdgeEventData.eventTime = time;
				}

				eventTimeListener.UpdateEdgeEvent(prevVertexData.prevEdgeIndex);
				eventTimeListener.UpdateEdgeEvent(nextVertexData.nextEdgeIndex);
			}

			// now update velocity for all adjacent vertices, but only if we didn't finish with a nook
			else if (updateEdgeEvents)
			{
				UpdateVertexData(newVertexIndex, time);
				UpdateVertexData(prevEdge.prevVertexIndex, time);
				UpdateVertexData(nextEdge.nextVertexIndex, time);

				// otherwise, we continue the process of updating event times for the edges
				UpdateEventTime(vertexDatas[newVertexIndex].prevEdgeIndex, time);
				eventTimeListener.UpdateEdgeEvent(vertexDatas[newVertexIndex].prevEdgeIndex);
				UpdateEventTime(vertexDatas[newVertexIndex].nextEdgeIndex, time);
				eventTimeListener.UpdateEdgeEvent(vertexDatas[newVertexIndex].nextEdgeIndex);
			}

			return newVertexIndex;
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

			// first, the easy part - calculate the edge event
			UpdateEdgeEventTime(ref edge, ref eventData, time);

			// then we see if this edge is being split by a vertex
			UpdateSplitEventTime(edgeIndex, ref eventData, time);
		}

		/*private void UpdateEdgeEventTime(ref Edge edge, ref EdgeEvent eventData, float currentTime)
		{
			// we store the edge event time in the prevVertex associated with the edge!
			float2 prevVertex = GetVertexPosAtTime(edge.prevVertexIndex, currentTime);
			float2 nextVertex = GetVertexPosAtTime(edge.nextVertexIndex, currentTime);
			ref VertexData prevData = ref vertexDatas[edge.prevVertexIndex];
			ref VertexData nextData = ref vertexDatas[edge.nextVertexIndex];

			// if the velocity is parallel to the line, we don't have an edge event
			float2 edgeDir = nextVertex - prevVertex;
			if (Geometry.Angle(edgeDir, math.normalize(prevData.velocity)) < Geometry.EPS) return;
			if (Geometry.Angle(edgeDir, math.normalize(nextData.velocity)) < Geometry.EPS) return;
			if (math.lengthsq(prevData.velocity) < Geometry.EPSSQ) return;
			if (math.lengthsq(nextData.velocity) < Geometry.EPSSQ) return;

			if (Geometry.GetLineIntersection(prevVertex, prevVertex + prevData.velocity, nextVertex, nextVertex + nextData.velocity, out float t0, out float t1))
			{
				if (t0 > -Geometry.EPS && t1 > -Geometry.EPS)
				{
					eventData.eventType = EventType.Edge;
					eventData.eventPos = prevVertex + prevData.velocity * t0;

					float2 projPoint = Geometry.ProjectPointOnLine(prevVertex, nextVertex, eventData.eventPos, out float t);
					float distance = math.distance(projPoint, eventData.eventPos);
					eventData.eventTime = currentTime + distance;
				}
			}
		}*/

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
					// the REMAINING time it takes to get to the collapse point.
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
			if (time < vertexData.creationTime - Geometry.EPS) throw new ArgumentException($"Time {time} lies in the past of the creation point of vertex {vertex} ({vertexData.creationTime}), this should not be possible.");
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

			// we can now determine the position where the edge meets the reflex vertex
			float distanceAlongVertexLine = L * V;
			float2 eventPos = reflexVertex + math.normalize(reflexVertex) * distanceAlongVertexLine;

			// because lines always move at constant speed, the time is identical to the move distance of the edge
			float time = L;
			float eventTime = currentTime + time;

			// if the reflex vertex collides with us after our collapse, we don't bother
			if (edgeEvent.eventType == EventType.Edge && eventTime > edgeEvent.eventTime - Geometry.EPS) return;

			// we have multiple splits at the same time - we just skip this one and do them one by one
			if (edgeEvent.eventType == EventType.Split && math.abs(eventTime - edgeEvent.eventTime) < Geometry.EPS) return;

			// if at this point we still have an event happening at the same time, we don't support it
			if (math.abs(eventTime - edgeEvent.eventTime) < Geometry.EPS) throw new NotSupportedException();

			// now finally, we need to make sure that this still falls into the edge once it moved that far
			float2 prevVertexAtEventTime = GetVertexPosAtTime(edge.prevVertexIndex, eventTime);
			float2 nextVertexAtEventTime = GetVertexPosAtTime(edge.nextVertexIndex, eventTime);
			Geometry.ProjectPointOnLine(prevVertexAtEventTime, nextVertexAtEventTime, eventPos, out float tLine);
			if (tLine < 0 || tLine > 1) return;

			// there was another event earlier
			if (edgeEvent.eventTime < eventTime) return;

			// we have a split event!
			edgeEvent.eventType = EventType.Split;
			edgeEvent.eventTime = eventTime;
			edgeEvent.eventPos = eventPos;
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
			Debug.Log("Vertex #" + vertexIndex + " has dir " + dir + " and speed " + speed + " based on " + moveDir1 + " and " + moveDir2 + " and final dir has angle " + (Geometry.GetAngle(dir) * Mathf.Rad2Deg) + ", angle between two adjacent dirs is " + (Geometry.SignedAngle(prevToCurr, nextToCurr) * Mathf.Rad2Deg));
			return velocity;
		}

		public float GetVertexTime(int vertexIndex)
		{
			return vertexDatas[vertexIndex].creationTime;
		}

		public void DrawGizmos(float time)
		{
			Gizmos.color = Color.cyan;

			for (int i = 0; i < nEdges; ++i)
			{
				ref var edge = ref edges[i];
				if (vertexDatas[edge.prevVertexIndex].inWavefront && vertexDatas[edge.nextVertexIndex].inWavefront)
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
				var vertex = vertices[i];
				ref var vertexData = ref vertexDatas[i];
				if (!vertexData.inWavefront) continue;

				var endPoint = vertex + math.normalize(vertexData.velocity) * 0.3f;

				Gizmos.DrawLine(new Vector3(vertex.x, 0.01f, vertex.y), new Vector3(endPoint.x, 0.01f, endPoint.y));
			}

		}
	}
}