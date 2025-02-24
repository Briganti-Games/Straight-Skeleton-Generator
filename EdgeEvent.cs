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
	public struct EdgeEvent
	{
		public int queueId;

		public EventType eventType;
		public float eventTime;
		public float2 eventPos;

		public int firstSplitReflexVertexIndex;


		public void Reset()
		{
			eventType = EventType.None;
			eventTime = float.MaxValue;
			firstSplitReflexVertexIndex = -1;
		}

		public override string ToString()
		{
			if (eventType == EventType.NotInWavefront) return "";
			return $"{eventType} event at time {eventTime} and pos {eventPos}";
		}
	}
}