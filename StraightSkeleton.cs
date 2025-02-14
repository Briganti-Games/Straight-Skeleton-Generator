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
		public readonly float[] vertexTimes;

		public StraightSkeleton(int maxVertices, int maxEdges) : base(maxVertices, maxEdges)
		{
			vertexTimes = new float[maxVertices];
		}

		public int AddVertex(float2 pos, float time)
		{
			int vertexIndex = base.AddVertex(pos);
			vertexTimes[vertexIndex] = time;
			return vertexIndex;
		}
	}
}