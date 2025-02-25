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
	public enum SplitPoint
	{
		PrevVertex,
		Edge,
		NextVertex,
	}

	public struct EdgeEvent
	{
		public int queueId;

		public EventType eventType;
		public float eventTime;
		public float2 eventPos;

		public SplitPoint splitPoint;
		public int reflexVertexIndex;

		public void Reset()
		{
			eventType = EventType.None;
			eventTime = float.MaxValue;
			reflexVertexIndex = -1;
		}

		public override string ToString()
		{
			if (eventType == EventType.NotInWavefront) return "";
			return $"{eventType} event at time {eventTime} and pos {eventPos}";
		}
	}
}