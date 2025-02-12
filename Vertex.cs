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
		public float2 pos => key.pos;
		public readonly float2 velocity;
		public WavefrontVertexType type;
		public readonly List<Edge> adjacentEdges;

		public readonly WavefrontVertexKey key;



		public Vertex(float2 pos, float2 velocity, WavefrontVertexType type)
		{
			key = new WavefrontVertexKey(pos);
			this.velocity = velocity;
			this.type = type;
			this.adjacentEdges = new List<Edge>();
		}

		public bool IsMoving()
		{
			return type != WavefrontVertexType.SteinerResting && type != WavefrontVertexType.SteinerMulti;
		}

		public bool IsWavefrontVertex()
		{
			return type == WavefrontVertexType.Convex || type == WavefrontVertexType.ConvexAndSteiner || type == WavefrontVertexType.Reflex || type == WavefrontVertexType.SteinerMoving;
		}

		public bool HasVelocity()
		{
			return math.lengthsq(velocity) > Geometry.EPS * Geometry.EPS;
		}

		public override string ToString()
		{
			return $"Vertex {type} at {pos} moving at velocity {velocity}";
		}
	}
}