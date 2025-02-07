using Briganti.StraightSkeletons.Priority_Queue;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Profiling;
using UnityEngine.Rendering;
using UnityEngine.UI;

namespace Briganti.StraightSkeletons
{
	public struct WavefrontEdgeSplitPoint : IComparable<WavefrontEdgeSplitPoint>
	{
		public int index;
		public float t;

		public WavefrontEdgeSplitPoint(int index, float t)
		{
			this.index = index;
			this.t = t;
		}

		public int CompareTo(WavefrontEdgeSplitPoint other)
		{
			return t.CompareTo(other.t);
		}
	}

	public class WavefrontEdge : FastPriorityQueueNode
	{
		public readonly int prevVertexIndex;
		public readonly int nextVertexIndex;

		public List<WavefrontEdgeSplitPoint> splitPoints;

		public float collapseTime { get; private set; }
		public float2 collapsePlace { get; private set; }

		public WavefrontEdge(int prevVertexIndex, int nextVertexIndex)
		{
			this.prevVertexIndex = prevVertexIndex;
			this.nextVertexIndex = nextVertexIndex;
		}

		public void UpdateCollapse(float collapseTime, float2 collapsePos)
		{
			this.collapseTime = collapseTime;
			this.collapsePlace = collapsePlace;
		}

		public override string ToString()
		{
			return $"Edge from {prevVertexIndex} to {nextVertexIndex} collapsing at {collapsePlace} at time {collapseTime}";
		}
	}
}