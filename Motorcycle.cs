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
	public struct Motorcycle
	{
		public float startTime;
		public readonly float2 startPoint;
		public readonly float2 velocity;

		public float crashTime;

		public Motorcycle(float startTime, float2 startPoint, float2 velocity)
		{
			this.startTime = startTime;
			this.startPoint = startPoint;
			this.velocity = velocity;

			this.crashTime = float.MaxValue;
		}

		public float2 GetPointOnTrace(float time)
		{
			if (time < startTime) throw new ArgumentException($"The motorcycle starting at {startTime} at {startPoint} with velocity {velocity} hasn't started yet at time {time}.");
			return startPoint + velocity * (time - startTime);
		}
	}
}