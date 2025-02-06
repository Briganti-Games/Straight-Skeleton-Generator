
using System;
using System.Collections.Generic;
using System.Linq;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

namespace Briganti.StraightSkeletons
{
	public class MotorcycleGraphGenerator
	{
		private readonly PolygonWavefront polygonWavefront;

		private Motorcycle[] motorcycles;
		private MotorcycleCrash[] motorcycleCrashes;
		private float maxTime = 0.0f;

		public MotorcycleGraphGenerator(PolygonWavefront polygonWavefront)
		{
			this.polygonWavefront = polygonWavefront;
		}

		public Motorcycle[] CalculateGraph()
		{

			// first, we generate motorcycles from all reflex angles
			motorcycles = GenerateMotorcycles();

			// now collide the motorcycles with all edges of the polygon - this is always the final stop of the motorcycle
			CollideMotorcyclesWithPolygon();

			// we now compute all the collisions between motorcycles
			CollideWithOtherMotorcycles();

			// all done!
			return motorcycles;
		}

		private Motorcycle[] GenerateMotorcycles()
		{
			List<Motorcycle> motorcycles = new List<Motorcycle>();

			// go over all contours and generate the motorcycles
			GenerateMotorcyclesForContour(polygonWavefront.outerWavefront, motorcycles);
			for (int i = 0; i < polygonWavefront.innerWavefronts.Count; ++i)
			{
				var innerWavefront = polygonWavefront.innerWavefronts[i];
				GenerateMotorcyclesForContour(innerWavefront, motorcycles);
			}

			return motorcycles.ToArray();
		}

		private void GenerateMotorcyclesForContour(Wavefront wavefront, List<Motorcycle> motorcycles)
		{
			var vertices = wavefront.vertices;
			for (int i = 0; i < vertices.Length; ++i)
			{
				//WavefrontVertex prev = vertices[(i - 1 + vertices.Length) % vertices.Length];
				WavefrontVertex curr = vertices[i];
				//WavefrontVertex next = vertices[(i + 1) % vertices.Length];

				if (curr.type == WavefrontVertexType.Reflex)
				{
					Motorcycle motorcycle = new Motorcycle(0.0f, curr.pos, curr.velocity);
					motorcycles.Add(motorcycle);
				}
			}
		}

		private void CollideMotorcyclesWithPolygon()
		{
			for (int i = 0; i < motorcycles.Length; ++i)
			{
				CollideMotorcycleWithPolygon(ref motorcycles[i]);
			}
		}

		private void CollideMotorcycleWithPolygon(ref Motorcycle motorcycle)
		{
			CollideMotorcycleWithWavefront(ref motorcycle, polygonWavefront.outerWavefront.vertices);
			for (int i = 0; i < polygonWavefront.innerWavefronts.Count; ++i)
			{
				CollideMotorcycleWithWavefront(ref motorcycle, polygonWavefront.innerWavefronts[i].vertices);
			}
		}

		private void CollideMotorcycleWithWavefront(ref Motorcycle motorcycle, WavefrontVertex[] vertices)
		{
			float2 motorcyclePoint0 = motorcycle.startPoint;
			float2 motorCyclePoint1 = motorcycle.startPoint + motorcycle.velocity;

			for (int i = 0; i < vertices.Length; ++i)
			{
				float2 edgePoint0 = vertices[i].pos;
				float2 edgePoint1 = vertices[(i + 1) % vertices.Length].pos;

				if (Geometry.GetLineIntersection(motorcyclePoint0, motorCyclePoint1, edgePoint0, edgePoint1, out float motorcycleTime, out float edgeTime))
				{
					// we hit the actual edge
					if (0.0f <= edgeTime && edgeTime <= 1.0f)
					{
						// don't hit adjacent edges that we originate from - so we make sure that we have at least moved a LITTLE forward
						if (motorcycleTime > 0.001f)
						{
							motorcycleTime += motorcycle.startTime;
							motorcycle.crashTime = Mathf.Min(motorcycle.crashTime, motorcycleTime);
							maxTime = Mathf.Max(maxTime, motorcycle.crashTime);
						}
					}
				}
			}
		}

		private void CollideWithOtherMotorcycles()
		{
			// we create a array of crashes per motorcycle
			List<MotorcycleCrash> crashes = new List<MotorcycleCrash>();
			List<MotorcycleCrash> realCrashes = new List<MotorcycleCrash>();
			List<Motorcycle> newMotorcycles = new List<Motorcycle>();
			int nMotorcycles = motorcycles.Length;
			int nUnprocessedMotorcycles = nMotorcycles;

			while (nUnprocessedMotorcycles > 0)
			{
				// flip the two lists
				/*var newRealCrashes = crashes;
				crashes = realCrashes;
				realCrashes = newRealCrashes;
				realCrashes.Clear();*/
				crashes.Clear();
				realCrashes.Clear();

				// we move from back to front in the list because the new motorcycles are appended at the back!
				for (int i = 0; i < nUnprocessedMotorcycles; ++i)
				{
					int motorcycle0Index = nMotorcycles - 1 - i;
					for (int motorcycle1Index = 0; motorcycle1Index < motorcycle0Index; ++motorcycle1Index)
					{
						ref Motorcycle motorcycle0 = ref motorcycles[motorcycle0Index];
						ref Motorcycle motorcycle1 = ref motorcycles[motorcycle1Index];
						float2 p00 = motorcycle0.startPoint;
						float2 p01 = p00 + motorcycle0.velocity;
						float2 p10 = motorcycle1.startPoint;
						float2 p11 = p10 + motorcycle1.velocity;

						if (Geometry.GetLineIntersection(p00, p01, p10, p11, out float t0, out float t1))
						{
							float time0 = motorcycle0.startTime + t0;
							float time1 = motorcycle1.startTime + t1;

							// we only move forward
							if (motorcycle0.startTime < time0 && time0 <= motorcycle0.crashTime && motorcycle1.startTime < time1 && time1 <= motorcycle1.crashTime)
							{
								// we crash at the exact same time - both crashes are added
								if (Mathf.Abs(time0 - time1) < Geometry.EPS)
								{
									crashes.Add(new MotorcycleCrash(motorcycle0Index, motorcycle1Index, time0, time1, motorcycle0.GetPointOnTrace(time0)));
									crashes.Add(new MotorcycleCrash(motorcycle1Index, motorcycle0Index, time1, time0, motorcycle1.GetPointOnTrace(time1)));
								}

								// when they meet, motorcycle0 took less time to get there, so it mans motorcycle1 crashed against its trace
								else if (time0 < time1) crashes.Add(new MotorcycleCrash(motorcycle1Index, motorcycle0Index, time1, time0, motorcycle1.GetPointOnTrace(time1)));

								// motorcycle1 took less time to get here, so has passed, and motorcycle0 crashed against its trace
								else crashes.Add(new MotorcycleCrash(motorcycle0Index, motorcycle1Index, time0, time1, motorcycle0.GetPointOnTrace(time0)));
							}
						}
					}
				}

				// now we sort the crashes by time of when they occur
				crashes.Sort();
				int nCrashes = crashes.Count;

				// we now prune any invalid crashes by only counting the first crash that occurs for each motorcycle
				for (int i = 0; i < nCrashes; ++i)
				{
					MotorcycleCrash crash = crashes[i];

					// this motorcycle previously crashed - skip this crash
					if (!crash.valid) continue;

					// disable any later crashes of this motorcycle
					//DisableLaterCrashes(crash.motorcycleIndex, i + 1, crash.time);

					// we go over all other crashes for this motorcycle and we set them invalid - they don't matter anymore since we found a real crash
					for (int j = i + 1; j < nCrashes; ++j)
					{
						MotorcycleCrash laterCrash = crashes[j];

						// this is a crash of the same motorcycle later down its lifecycle - this crash will never occur because it already crashed
						if (laterCrash.motorcycleIndex == crash.motorcycleIndex) laterCrash.valid = false;

						// another motorcycle crashed against our trace LATER than now - we also delete it
						// we make sure that the crash is LATER because it might occur at the exact same time, in which case we want the crash to be processed later
						else if (laterCrash.motorcycleTraceIndex == crash.motorcycleIndex && laterCrash.traceTime > crash.time + Geometry.EPS)
						{
							laterCrash.valid = false;
						}

						crashes[j] = laterCrash;
					}

					ref Motorcycle motorcycle = ref motorcycles[crash.motorcycleIndex];
					motorcycle.crashTime = Mathf.Min(motorcycle.crashTime, crash.time);

					realCrashes.Add(crash);
				}

				// now we go over all real crashes, and we detect simultaneous crash events - in this case, we might need to spawn new motorcycles
				nUnprocessedMotorcycles = 0;
				newMotorcycles.Clear();
				int nSimultaneousCrashes = 1;
				int nRealCrashes = realCrashes.Count;
				for (int i = 0; i < nRealCrashes - 1; ++i)
				{

					// crash at the same time & place
					if (realCrashes[i].CompareTo(realCrashes[i + 1]) == 0)
					{
						++nSimultaneousCrashes;
					}

					// we have a crash that is NOT at the same time & place - in that case, we see if we need to launch additional motorcycles
					else
					{
						if (nSimultaneousCrashes > 1)
						{
							nUnprocessedMotorcycles += AddNewMotorcycles(realCrashes, i - nSimultaneousCrashes + 1, nSimultaneousCrashes, newMotorcycles);
							nSimultaneousCrashes = 1;
						}
					}
				}

				// we might still have a simultaneous crash going at the end of the line
				if (nSimultaneousCrashes > 1)
				{
					nUnprocessedMotorcycles += AddNewMotorcycles(realCrashes, realCrashes.Count - 1 - nSimultaneousCrashes + 1, nSimultaneousCrashes, newMotorcycles);
				}

				// append the new motorcycles to the old ones
				Motorcycle[] allMororcycles = new Motorcycle[nMotorcycles + nUnprocessedMotorcycles];
				Array.Copy(motorcycles, allMororcycles, nMotorcycles);
				for (int i = 0; i < newMotorcycles.Count; ++i)
				{
					allMororcycles[nMotorcycles + i] = newMotorcycles[i];

					// make sure we also calculate wall collisions for these new motorcycles!
					ref Motorcycle motorcycle = ref allMororcycles[nMotorcycles + i];
					CollideMotorcycleWithPolygon(ref motorcycle);
				}

				motorcycles = allMororcycles;
				nMotorcycles = motorcycles.Length;
			}
		}


		private List<Vector2> incomingMotorcycleDirs = new List<Vector2>();
		private int AddNewMotorcycles(List<MotorcycleCrash> crashes, int startIndex, int nSimultaneousCrashes, List<Motorcycle> newMotorcycles)
		{
			int nNewMotorcycles = 0;
			float2 p = crashes[startIndex].crashPos;
			float crashTime = crashes[startIndex].time;

			// collect and sort all incoming motorcycles
			incomingMotorcycleDirs.Clear();
			for (int i = startIndex; i < startIndex + nSimultaneousCrashes; ++i)
			{
				var motorcycle = motorcycles[crashes[i].motorcycleIndex];
				incomingMotorcycleDirs.Add(motorcycle.velocity);
			}

			// sort them by angle, counter-clockwise
			incomingMotorcycleDirs.Sort((dir1, dir2) => Geometry.GetAngle(dir1).CompareTo(Geometry.GetAngle(dir2)));

			// now go over each couple, and see if there's a reflex angle between them - if there is, we launch a new motorcycle in the middle between them
			for (int i = 0; i < incomingMotorcycleDirs.Count; ++i) {
				float2 dir0 = incomingMotorcycleDirs[i];
				float2 dir1 = incomingMotorcycleDirs[(i + 1) % incomingMotorcycleDirs.Count];

				float2 prev = p - dir1;
				float2 curr = p;
				float2 next = p - dir0;

				// we flip them because we need to find the relfex 
				if (Geometry.IsRelfexVertex(prev, curr, next))
				{
					float2 newVelocity = Wavefront.CalculateVelocity(prev, curr, next);
					var newMotorcycle = new Motorcycle(crashTime, p, newVelocity);
					newMotorcycles.Add(newMotorcycle);
					++nNewMotorcycles;
				}
			}

			return nNewMotorcycles;

		}
	}
}