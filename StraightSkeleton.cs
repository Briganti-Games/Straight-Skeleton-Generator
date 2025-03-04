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
	public class StraightSkeleton : VertexGraph
	{
		public float[] vertexTimes { get; private set; }
		public int[] wavefrontToStraightSkeletonVertexIndices { get; private set; }

		public float maxTime { get; private set; } = float.MinValue;

		public StraightSkeleton(int maxVertices, int maxEdges) : base(maxVertices, maxEdges)
		{
			vertexTimes = new float[maxVertices];
		}

		public int AddVertex(float2 pos, float time)
		{
			int vertexIndex = base.AddVertex(pos);
			vertexTimes[vertexIndex] = time;
			maxTime = Mathf.Max(maxTime, time);
			return vertexIndex;
		}

		protected override void IncreaseVertexCapacity(int nVertices)
		{
			base.IncreaseVertexCapacity(nVertices);

			int oldLength = vertexTimes.Length;

			float[] newVertexTimes = new float[nVertices];
			Array.Copy(vertexTimes, newVertexTimes, oldLength);
			vertexTimes = newVertexTimes;
		}
		/*
		public List<List[]> GetIndices()
		{
			Profiler.BeginSample("Biganti.StraightSkeleton.GetIndices");
			List<List<int>[]> polygons = new List<List<int>[]>();

			bool[] visited = new bool[nEdges];
			for (int edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex)
			{
				if (visited[edgeIndex]) continue;

				visited[edgeIndex] = true;
				ref var edge = ref edges[edgeIndex];

				List<int> polygon = new List<int>();
				polygon.Add(edge.prevVertexIndex);
				polygon.Add(edge.nextVertexIndex);

				ref var nextEdge = ref edges[ver

				
			}


			Profiler.EndSample();

			return polygons;
		}*/
	}
}