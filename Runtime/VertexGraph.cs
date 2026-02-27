using System;
using Unity.Mathematics;
using UnityEngine;

namespace Briganti.StraightSkeletonGeneration
{
	public class VertexGraph
	{
		public float2[] vertices { get; private set; }
		public int nVertices { get; private set; } = 0;
		public int maxVertices { get; private set; }

		public Edge[] edges { get; private set; }
		public int nEdges { get; private set; } = 0;
		public int maxEdges { get; private set; }


		public VertexGraph()
		{
			// we initialize later
		}

		public VertexGraph(int maxVertices, int maxEdges)
		{
			Initialize(maxVertices, maxEdges);
		}

		protected void Initialize(int maxVertices, int maxEdges)
		{
			this.maxVertices = maxVertices;
			this.maxEdges = maxEdges;

			// this can be improved and the upper limit can be calculated EXACTLY - see the lemmas in Stefan Huber's PhD
			vertices = new float2[maxVertices];
			edges = new Edge[maxEdges];
		}

		public int AddVertex(float2 pos)
		{
			if (nVertices == maxVertices) IncreaseVertexCapacity(maxVertices * 2);

			// we need to calculate the convex/reflex state later when the adjacent edges are known!!
			vertices[nVertices++] = new float2(pos);
			return nVertices - 1;
		}

		public int AddEdge(int prevVertexIndex, int nextVertexIndex)
		{
			return AddEdge(new Edge(prevVertexIndex, nextVertexIndex));
		}

		public int AddEdge(Edge edge)
		{
			if (nEdges == maxEdges) IncreaseEdgeCapacity(maxEdges * 2);

			edges[nEdges++] = edge;
			return nEdges - 1;
		}

		protected virtual void IncreaseVertexCapacity(int nVertices)
		{
			int oldMaxVertices = maxVertices;
			maxVertices = nVertices;
			float2[] newVertices = new float2[maxVertices];
			Array.Copy(vertices, newVertices, oldMaxVertices);
			vertices = newVertices;
		}

		protected virtual void IncreaseEdgeCapacity(int nEdges)
		{
			int oldMaxEdges = maxEdges;
			maxEdges = nEdges;
			Edge[] newEdges = new Edge[maxEdges];
			Array.Copy(edges, newEdges, oldMaxEdges);
			edges = newEdges;
		}
	}
}