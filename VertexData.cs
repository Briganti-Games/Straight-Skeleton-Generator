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
	public struct VertexData
	{
		public float2 velocity;
		public WavefrontVertexType type;
		public float creationTime;

		public int prevVertexIndex;
		public int nextVertexIndex;

		public int prevEdgeIndex;
		public int nextEdgeIndex;

		public bool inWavefront;

		// this is set by the reflex vertex itself
		public SplitPoint splitPoint;
		public bool partOfSplitEvent;
		public float splitTime;
		public int splitEdge;

		// this is set after assigning the vertex to an edge event
		//public int nextSplitReflexVertexIndex;
		

		public void UpdateConnections(int prevVertexIndex, int nextVertexIndex, int prevEdgeIndex, int nextEdgeIndex)
		{
			this.prevVertexIndex = prevVertexIndex;
			this.nextVertexIndex = nextVertexIndex;
			this.prevEdgeIndex = prevEdgeIndex;
			this.nextEdgeIndex = nextEdgeIndex;
			this.inWavefront = true;
		}

		public override string ToString()
		{
			if (!inWavefront) return "";
			return $"{type} vertex with velocity {velocity} created at {creationTime}, prev vertex {prevVertexIndex}, next vertex {nextVertexIndex}, still active {inWavefront}";
		}
	}
}
