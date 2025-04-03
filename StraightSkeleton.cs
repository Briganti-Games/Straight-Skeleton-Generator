using System;
using Unity.Mathematics;
using UnityEngine;

namespace Briganti.StraightSkeletonGeneration
{
	public class StraightSkeleton : VertexGraph
	{
		public float[] vertexTimes { get; private set; }

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
	}
}