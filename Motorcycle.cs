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
	public enum MotorcycleCrashType
	{
		Wall,
		Motorcycle,
		Escaped,
	}

	public struct Motorcycle
	{
		public int vertexIndex;

		public float startTime;
		public readonly float2 startPoint;
		public readonly float2 velocity;

		public float crashTime { get; private set; }
		public MotorcycleCrashType crashType { get; private set; }
		public int crashTargetIndex { get; private set; }

		public Motorcycle(int vertexIndex, float startTime, float2 startPoint, float2 velocity)
		{
			this.vertexIndex = vertexIndex;
			this.startTime = startTime;
			this.startPoint = startPoint;
			this.velocity = velocity;

			this.crashType = MotorcycleCrashType.Escaped;
			this.crashTime = float.MaxValue;
			this.crashTargetIndex = 0;
		}

		public void UpdateCrash(MotorcycleCrashType crashType, int crashTargetIndex, float crashTime)
		{
			this.crashType = crashType;
			this.crashTargetIndex = crashTargetIndex;
			this.crashTime = crashTime;
		}

		public float2 GetPointOnTrace(float time)
		{
			if (time < startTime) throw new ArgumentException($"The motorcycle starting at {startTime} at {startPoint} with velocity {velocity} hasn't started yet at time {time}.");
			return startPoint + velocity * (time - startTime);
		}

		public float2 getCrashPos()
		{
			return GetPointOnTrace(crashTime);
		}

		public override string ToString()
		{
			return $"Motorcycle starting at {startPoint} with velocity {velocity}, crashing against {crashType} at time {crashTime}";
		}
	}
}