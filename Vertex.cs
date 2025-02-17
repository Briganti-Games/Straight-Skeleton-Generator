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
		Nook,
	}

	public static class EventTypeExtensions
	{
		public static bool IsBatchEvent(this EventType eventType)
		{
			return eventType == EventType.Edge;
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