using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Profiling;
using UnityEngine.Rendering;
using UnityEngine.UI;

namespace Briganti.StraightSkeletonGeneration
{
	public class WavefrontGraph : VertexGraph
	{
		public readonly EdgeEvent[] edgeEvents;
		public readonly VertexData[] vertexDatas;


		public WavefrontGraph(PolygonWithHoles polygonWithHoles)
		{
			int nVertices = polygonWithHoles.outerContour.Length;
			for (int i = 0; i < polygonWithHoles.innerContours.Count; ++i)
			{
				nVertices += polygonWithHoles.innerContours[i].Length;
			}

			int nEdges = nVertices;

			// we now add all the edges that will be added for the straight skeleton - the exact number is known
			nVertices += nVertices - 2;
			int nArcs = ((2 * nVertices) - 3) * 2; // multiply by 2 because we will be duplicating both sides of an edge for the polygon on each side
			nEdges += nArcs;

			// this can be improved and the upper limit can be calculated EXACTLY - see the lemmas in Stefan Huber's PhD
			Initialize(nVertices, nEdges);

			// we store additional data for each edge and vertex
			edgeEvents = new EdgeEvent[nEdges];
			vertexDatas = new VertexData[nVertices];

			// add all contours
			AddContour(polygonWithHoles.outerContour);
			for (int i = 0; i < polygonWithHoles.innerContours.Count; ++i)
			{
				AddContour(polygonWithHoles.innerContours[i]);
			}

			// update all events
			for (int i = 0; i < nEdges; ++i)
			{
				UpdateEventTime(i);
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
				int vertexIndex = AddVertex(curr);

				int prevVertexIndex = vertexOffset + (i - 1 + contour.Length) % contour.Length;
				int nextVertexIndex = vertexOffset + (i + 1) % contour.Length;

				// set the vertex data
				int prevEdgeIndex = edgeOffset + (i - 1 + contour.Length) % contour.Length;
				int nextEdgeIndex = edgeOffset + (i + 1) % contour.Length;
				vertexDatas[i] = new VertexData(prevVertexIndex, nextVertexIndex, prevEdgeIndex, nextEdgeIndex);

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

		public void AddVertexToWavefrontAndRemoveEdge(float2 pos, float depth, int disappearedEdgeIndex)
		{
			// firstly, just add the vertex
			int newVertexIndex = AddVertex(pos);

			// now we cut out the disappeared edge from the graph and immediately connect the adjacent vertices
			ref Edge disappearedEdge = ref edges[disappearedEdgeIndex];
			ref VertexData prevVertexData = ref vertexDatas[disappearedEdge.prevVertexIndex];
			ref VertexData nextVertexData = ref vertexDatas[disappearedEdge.nextVertexIndex];

			if (prevVertexData.nextEdgeIndex != disappearedEdgeIndex) throw new ArgumentException($"Vertex {disappearedEdge.prevVertexIndex} is not properly linked to edge {disappearedEdgeIndex}");
			if (nextVertexData.prevEdgeIndex!= disappearedEdgeIndex) throw new ArgumentException($"Vertex {disappearedEdge.nextVertexIndex} is not properly linked to edge {disappearedEdgeIndex}");

			int prevEdgeIndex = prevVertexData.prevEdgeIndex;
			int nextEdgeIndex = nextVertexData.nextEdgeIndex;

			int prevVertexIndex = prevVertexData.prevVertexIndex;
			int nextVertexIndex = nextVertexData.nextVertexIndex;

			vertexDatas[newVertexIndex] = new VertexData(prevVertexIndex, nextVertexIndex, prevEdgeIndex, nextEdgeIndex);

			// now update velocity etc
			UpdateVertexData(newVertexIndex, depth);
		}

		private void UpdateVertexData(int vertexIndex, float depth) {

			ref Vertex vertex = ref vertices[vertexIndex];
			ref VertexData vertexData = ref vertexDatas[vertexIndex];
			ref Vertex prevVertex = ref vertices[vertexData.prevVertexIndex];
			ref Vertex nextVertex = ref vertices[vertexData.nextVertexIndex];

			float2 velocity = CalculateVelocity(prevVertex.pos, vertex.pos, nextVertex.pos);
			bool isReflex = Geometry.IsRelfexVertex(prevVertex.pos, vertex.pos, nextVertex.pos);
			var type = (isReflex ? WavefrontVertexType.Reflex : WavefrontVertexType.Convex);

			vertexData.velocity = velocity;
			vertexData.type = type;
			vertexData.depth = depth;
		}

		public void UpdateEventTime(int edgeIndex)
		{
			// default to max
			ref Edge edge = ref edges[edgeIndex];
			ref EdgeEvent eventData = ref edgeEvents[edgeIndex];
			eventData.eventTime = float.MaxValue;

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
			ref Vertex prevVertex = ref vertices[edge.prevVertexIndex];
			ref Vertex nextVertex = ref vertices[edge.nextVertexIndex];
			ref VertexData prevData = ref vertexDatas[edge.prevVertexIndex];
			ref VertexData nextData = ref vertexDatas[edge.nextVertexIndex];

			if (Geometry.GetLineIntersection(prevVertex.pos, prevVertex.pos + prevData.velocity, nextVertex.pos, nextVertex.pos + nextData.velocity, out float t0, out float t1))
			{
				if (t0 > 0 && t1 > 0)
				{
					eventData.eventType = EventType.Edge;
					eventData.eventPos = prevVertex.pos + prevData.velocity * t0;
					Geometry.ProjectPointOnLine(prevVertex.pos, nextVertex.pos, eventData.eventPos, out float time);
					eventData.eventTime = time;
				}
			}
		}

		private void UpdateSplitEventTime(int vertexIndex)
		{
			// TODO
		}

		public static float2 CalculateVelocity(in float2 prev, in float2 curr, in float2 next)
		{
			float2 prevToCurr = curr - prev;
			float2 nextToCurr = next - curr;

			// rotate them 90°
			float2 moveDir1 = math.normalize(Geometry.RotateMinus90Degrees(prevToCurr));
			float2 moveDir2 = math.normalize(Geometry.RotateMinus90Degrees(nextToCurr));

			float2 dir = math.normalize(moveDir1 + moveDir2);
			float speed = 2.0f / (math.dot(moveDir1, dir) + math.dot(moveDir2, dir));
			float2 velocity = dir * speed;

			return velocity;
		}
	}
}