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
	public enum EventType
	{
		None,
		Edge,
		Split,
	}

	public struct WavefrontVertexKey : IEquatable<WavefrontVertexKey>
	{
		public readonly float2 pos;

		public WavefrontVertexKey(float2 pos)
		{
			this.pos = pos;
		}

		public bool Equals(WavefrontVertexKey other)
		{
			int2 pos = Geometry.RoundToInt(this.pos / Geometry.EPS);
			int2 otherPos = Geometry.RoundToInt(other.pos / Geometry.EPS);
			return pos.x == otherPos.x && pos.y == otherPos.y;
		}

		public override bool Equals(object obj)
		{
			if (!(obj is Vertex)) return false;
			return ((Vertex)obj).Equals(this);
		}

		public override int GetHashCode()
		{
			return pos.GetHashCode();
		}
	}

	public struct Vertex
	{
		public float2 pos;

		public Vertex(float2 pos)
		{
			this.pos = pos;
		}

		public override string ToString()
		{
			return $"{pos}";
		}
	}
}