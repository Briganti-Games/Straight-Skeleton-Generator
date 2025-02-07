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
	public struct WavefrontVertex : IEquatable<WavefrontVertex>
	{
		public readonly float2 pos;
		public readonly float2 velocity;
		public WavefrontVertexType type;

		public WavefrontVertex(float2 pos, float2 velocity, WavefrontVertexType type)
		{
			this.pos = pos;
			this.velocity = velocity;
			this.type = type;
		}

		public bool Equals(WavefrontVertex other)
		{
			int2 pos = Geometry.RoundToInt(this.pos / Geometry.EPS);
			int2 otherPos = Geometry.RoundToInt(other.pos / Geometry.EPS);
			return pos.x == otherPos.x && pos.y == otherPos.y;
		}

		public override bool Equals(object obj)
		{
			if (!(obj is WavefrontVertex)) return false;
			return ((WavefrontVertex)obj).Equals(this);
		}

		public override int GetHashCode()
		{
			return pos.GetHashCode();
		}

		public override string ToString()
		{
			return $"Vertex {type} at {pos} moving at velocity {velocity}";
		}
	}
}