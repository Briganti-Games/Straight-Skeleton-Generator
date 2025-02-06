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
	public class WavefrontEdge : FastPriorityQueueNode
	{
		public readonly int prevVertexIndex;
		public readonly int nextVertexIndex;

		public float collapseTime { get; private set; }

		public WavefrontEdge(int prevVertexIndex, int nextVertexIndex)
		{
			this.prevVertexIndex = prevVertexIndex;
			this.nextVertexIndex = nextVertexIndex;
		}
	}
}