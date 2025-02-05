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
	public struct WavefrontVertex
	{
		public readonly float2 pos;
		public readonly float2 velocity;
		public readonly WavefrontVertexType type;

		public WavefrontVertex(float2 pos, float2 velocity, WavefrontVertexType type)
		{
			this.pos = pos;
			this.velocity = velocity;
			this.type = type;
		}
	}
}