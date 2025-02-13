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
		public float depth;

		public int prevVertexIndex;
		public int nextVertexIndex;

		public int prevEdgeIndex;
		public int nextEdgeIndex;

		public VertexData( int prevVertexIndex, int nextVertexIndex, int prevEdgeIndex, int nextEdgeIndex)
		{
			this.type = WavefrontVertexType.Unknown;
			this.velocity = float2.zero;
			this.depth = 0;
			this.prevVertexIndex = prevVertexIndex;
			this.nextVertexIndex = nextVertexIndex;
			this.prevEdgeIndex = prevEdgeIndex;
			this.nextEdgeIndex = nextEdgeIndex;
		}
	}
}