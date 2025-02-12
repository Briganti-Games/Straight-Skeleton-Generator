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

namespace Briganti.StraightSkeletons
{
	public class VertexGraph
	{
		public Vertex[] vertices { get; private set; }
		public int nVertices { get; private set; } = 0;
		public List<Edge> edges { get; private set; } = new List<Edge>();

		private Dictionary<WavefrontVertexKey, int> vertexIndices = new Dictionary<WavefrontVertexKey, int>();


		public VertexGraph(Vertex[] vertices, int nVertices, List<Edge> edges)
		{
			this.vertices = new Vertex[nVertices];
			this.nVertices = nVertices;
			for (int i = 0; i < nVertices; ++i)
			{
				this.vertices[i] = new Vertex(vertices[i].pos, vertices[i].velocity, vertices[i].type);
			}
			this.edges = edges;
		}

		public VertexGraph(PolygonWavefront polygonWavefront)
		{
			int vertexOffset = 0;
			int edgeOffset = 0;

			vertices = new Vertex[polygonWavefront.wavefronts.Sum(wavefront => wavefront.vertices.Length)];
			for (int wavefrontIndex = 0; wavefrontIndex < polygonWavefront.wavefronts.Count; ++wavefrontIndex)
			{
				var wavefront = polygonWavefront.wavefronts[wavefrontIndex];
				foreach (Vertex vertex in wavefront.vertices)
				{
					AddVertex(vertex);
				}

				for (int i = 0; i < wavefront.edges.Count; ++i)
				{
					var edge = wavefront.edges[i];
					AddEdge(new Edge(edge.prevVertexIndex + vertexOffset, edge.nextVertexIndex + vertexOffset));
				}

				vertexOffset += wavefront.vertices.Length;
				edgeOffset += wavefront.edges.Count;
			}
		}

		public VertexGraph(float2[] contour)
		{

			// convert all vertices in the contour to wavefront vertices
			vertices = new Vertex[contour.Length];
			edges = new List<Edge>(contour.Length);
			for (int i = 0; i < contour.Length; ++i)
			{
				float2 prev = contour[(i - 1 + contour.Length) % contour.Length];
				float2 curr = contour[i];
				float2 next = contour[(i + 1) % contour.Length];

				float2 velocity = CalculateVelocity(prev, curr, next);

				bool isReflex = Geometry.IsRelfexVertex(prev, curr, next);
				int prevEdgeIndex = (i - 1 + contour.Length) % contour.Length;
				int nextEdgeIndex = i;
				vertices[i] = new Vertex(curr, velocity, isReflex ? WavefrontVertexType.Reflex : WavefrontVertexType.Convex);
				edges.Add(new Edge(i, (i + 1) % contour.Length));
			}
		}

		public int AddVertex(in Vertex vertex)
		{
			if (vertices.Length == nVertices)
			{
				Vertex[] newVertices = new Vertex[nVertices * 2];
				Array.Copy(vertices, newVertices, nVertices);
				vertices = newVertices;
			}

			if (!vertexIndices.TryGetValue(vertex.key, out int vertexIndex))
			{
				vertexIndex = nVertices;
				vertices[nVertices++] = vertex;
				vertexIndices[vertex.key] = vertexIndex;
			}
			return vertexIndex;
		}

		public void AddEdge(Edge edge)
		{
			edges.Add(edge);

			// also make sure the vertices remember this connection
			ref Vertex prevVertex = ref vertices[edge.prevVertexIndex];
			prevVertex.adjacentEdges.Add(edge);

			ref Vertex nextVertex = ref vertices[edge.nextVertexIndex];
			nextVertex.adjacentEdges.Add(edge);
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

		public float2 GetMotorcycleDir(Edge edge, int vertexIndex)
		{
			ref var vertex = ref vertices[vertexIndex];
			if (vertex.type != WavefrontVertexType.SteinerMoving) throw new ArgumentException($"Vertex {vertex} is not a moving steiner vertex, so we can't calculate the motorcycle dir.");

			Edge motorcycleEdge;
			if (vertexIndex == edge.prevVertexIndex) {
				motorcycleEdge = edge.prevCounterClockwiseEdge;
			}
			else {
				motorcycleEdge = edge.nextClockwiseEdge;
			}

			return GetDirFrom(motorcycleEdge, vertexIndex);
		}

		public float GetAngleFrom(Edge edge, int vertexIndex)
		{
			float2 dir = GetDirFrom(edge, vertexIndex);
			return Geometry.GetAngle(dir);
		}

		public float2 GetDirFrom(Edge edge, int originVertexIndex)
		{
			ref Vertex prevVertex = ref vertices[edge.prevVertexIndex];
			ref Vertex nextVertex = ref vertices[edge.nextVertexIndex];

			if (edge.prevVertexIndex == originVertexIndex)
			{
				return math.normalize(nextVertex.pos - prevVertex.pos);
			}
			else
			{
				return math.normalize(prevVertex.pos - nextVertex.pos);
			}
		}

		public float2 GetDir(Edge edge) {
			return GetDirFrom(edge, edge.prevVertexIndex);
		}

		public float2 GetWavefrontDir(Edge edge)
		{
			if (!IsWavefrontEdge(edge)) throw new ArgumentException($"Edge {edge} is not a wavefront edge so it does not have a wavefront dir!");
			return Geometry.RotateMinus90Degrees(GetDir(edge));
		}

		public bool IsWavefrontEdge(Edge edge)
		{
			ref Vertex prevVertex = ref vertices[edge.prevVertexIndex];
			ref Vertex nextVertex = ref vertices[edge.nextVertexIndex];

			return prevVertex.IsWavefrontVertex() && nextVertex.IsWavefrontVertex();
		}
	}
}