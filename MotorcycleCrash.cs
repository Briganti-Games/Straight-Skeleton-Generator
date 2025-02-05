
using System;
using Unity.Mathematics;

namespace Briganti.StraightSkeletons
{
	public struct MotorcycleCrash : IComparable<MotorcycleCrash>
	{
		public int motorcycleIndex;
		public int motorcycleTraceIndex;
		public float time;
		public float traceTime;
		public bool valid;
		public float2 crashPos;

		public MotorcycleCrash(int motorcycleIndex, int motorcycleTraceIndex, float time, float traceTime, float2 crashPos)
		{
			this.motorcycleIndex = motorcycleIndex;
			this.motorcycleTraceIndex = motorcycleTraceIndex;
			this.time = time;
			this.traceTime = traceTime;
			this.crashPos = crashPos;
			this.valid = true;
		}

		public int CompareTo(MotorcycleCrash other)
		{
			// crashed at the same time
			if (math.abs(time - other.time) < Geometry.EPS)
			{
				// now sort by the crash position, so that consecutive crashes at the same position and time are sorted next to each other
				if (math.abs(crashPos.x - other.crashPos.x) < Geometry.EPS)
				{
					if (math.abs(crashPos.y - other.crashPos.y) < Geometry.EPS)
					{
						return 0;
					}
					else
					{
						return crashPos.y.CompareTo(other.crashPos.y);
					}
				}
				else
				{
					return crashPos.x.CompareTo(other.crashPos.x);
				}
			}
			else
			{
				return time.CompareTo(other.time);
			}
		}
	}
}