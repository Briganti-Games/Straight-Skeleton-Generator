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
		SplitAtVertex,
		Nook,
		NotInWavefront,
	}

	public static class EventTypeExtensions
	{
		public static bool IsBatchEvent(this EventType eventType)
		{
			return eventType == EventType.Edge || eventType == EventType.SplitAtVertex;
		}

		public static bool IsSplitEvent(this EventType eventType)
		{
			return eventType == EventType.Split || eventType == EventType.SplitAtVertex;
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