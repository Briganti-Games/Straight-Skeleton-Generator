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
				UpdateEventTime(i);
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

			// now update velocity for all adjacent vertices
			if (updateEdgeEvents)
			{
				UpdateVertexData(newVertexIndex, time);
				UpdateVertexData(prevEdge.prevVertexIndex, time);
				UpdateVertexData(nextEdge.nextVertexIndex, time);

				// otherwise, we continue the process of updating event times for the edges
				UpdateEventTime(vertexDatas[newVertexIndex].prevEdgeIndex);
				eventTimeListener.UpdateEdgeEvent(vertexDatas[newVertexIndex].prevEdgeIndex);
				UpdateEventTime(vertexDatas[newVertexIndex].nextEdgeIndex);
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

			vertexData.velocity = velocity;
			vertexData.type = type;
		}

		public void UpdateEventTime(int edgeIndex)
		{
			// default to max
			ref Edge edge = ref edges[edgeIndex];
			ref EdgeEvent eventData = ref edgeEvents[edgeIndex];

			eventData.Reset();

			// first, the easy part - calculate the edge event
			UpdateEdgeEventTime(ref edge, ref eventData);

			// if this vertex is reflex, we also update its split event time
			ref var vertexData = ref vertexDatas[edge.prevVertexIndex];
			if (vertexData.type == WavefrontVertexType.Unknown) throw new ArgumentException($"Vertex {edge.prevVertexIndex} was not properly initialized and added to the graph, because we don't know whether it is convex or reflex!");
			if (vertexData.type == WavefrontVertexType.Reflex)
			{
				UpdateSplitEventTime(edgeIndex);
			}
		}

		private void UpdateEdgeEventTime(ref Edge edge, ref EdgeEvent eventData)
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

		private void UpdateSplitEventTime(int vertexIndex)
		{
			// TODO
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
	}
}