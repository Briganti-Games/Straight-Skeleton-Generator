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

	public class Edge : FastPriorityQueueNode
	{
		public readonly int prevVertexIndex;
		public readonly int nextVertexIndex;

		public Edge prevClockwiseEdge, prevCounterClockwiseEdge;
		public Edge nextClockwiseEdge, nextCounterClockwiseEdge;

		public List<WavefrontEdgeSplitPoint> splitPoints;

		public float collapseTime { get; private set; }
		public float2 collapsePlace { get; private set; }

		public Edge(int prevVertexIndex, int nextVertexIndex)
		{
			this.prevVertexIndex = prevVertexIndex;
			this.nextVertexIndex = nextVertexIndex;
		}

		public Edge GetAdjacentEdge(int vertexIndex)
		{
			if (vertexIndex == prevVertexIndex) return prevClockwiseEdge;
			else if (vertexIndex == nextVertexIndex) return nextCounterClockwiseEdge;
			else throw new ArgumentException($"Vertex {vertexIndex} is not part of edge {this}.");
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