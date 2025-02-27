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
	public struct Edge
	{
		public int prevVertexIndex, nextVertexIndex;


		public Edge(int prevVertexIndex, int nextVertexIndex)
		{
			this.prevVertexIndex = prevVertexIndex;
			this.nextVertexIndex = nextVertexIndex;
		}

		public int GetPoint(SplitPoint point)
		{
			if (point == SplitPoint.PrevVertex) return prevVertexIndex;
			else if (point == SplitPoint.NextVertex) return nextVertexIndex;
			else throw new ArgumentException($"Split point {point} is not allowed.");
		}

		public override string ToString()
		{
			return $"Edge from vertices {prevVertexIndex} to {nextVertexIndex}";
		}
	}
}