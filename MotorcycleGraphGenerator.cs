
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
		private Motorcycle[] motorcycles;
		private MotorcycleCrash[] motorcycleCrashes;
		private float maxTime = 0.0f;

		private List<WavefrontVertex> vertices = new List<WavefrontVertex>();
		private List<WavefrontEdge> edges = new List<WavefrontEdge>();

		public MotorcycleGraphGenerator(PolygonWavefront polygonWavefront)
		{
			// first, we collect all existing vertices and edges from the wavefronts and throw them in a big bucket
			int vertexOffset = 0;
			int edgeOffset = 0;
			for (int wavefrontIndex = 0; wavefrontIndex < polygonWavefront.wavefronts.Count; ++wavefrontIndex)
			{
				var wavefront = polygonWavefront.wavefronts[wavefrontIndex];
				vertices.AddRange(wavefront.vertices);
				/*for (int i = 0; i < wavefront.vertices.Length; ++i)
				{
					var vertex = vertices[vertexOffset + i];
					vertex.prevEdgeIndex += edgeOffset;
					vertex.nextEdgeIndex += edgeOffset;
					vertices[vertexOffset + i] = vertex;
				}*/

				edges.AddRange(wavefront.edges);
				for (int i = 0; i < wavefront.edges.Length; ++i)
				{
					var edge = edges[edgeOffset + i];
					edges[vertexOffset + i] = new WavefrontEdge(edge.prevVertexIndex + vertexOffset, edge.nextVertexIndex + vertexOffset);
				}

				vertexOffset += wavefront.vertices.Length;
				edgeOffset += wavefront.edges.Length;
			}
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

			for (int i = 0; i < vertices.Count; ++i)
			{
				WavefrontVertex curr = vertices[i];

				if (curr.type == WavefrontVertexType.Reflex)
				{
					Motorcycle motorcycle = new Motorcycle(i, 0.0f, curr.pos, curr.velocity);
					motorcycles.Add(motorcycle);
				}
			}

			return motorcycles.ToArray();
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
			float2 motorcyclePoint0 = motorcycle.startPoint;
			float2 motorCyclePoint1 = motorcycle.startPoint + motorcycle.velocity;

			for (int i = 0; i < edges.Count; ++i)
			{
				var edge = edges[i];
				float2 edgePoint0 = vertices[edge.prevVertexIndex].pos;
				float2 edgePoint1 = vertices[edge.nextVertexIndex].pos;

				if (Geometry.GetLineIntersection(motorcyclePoint0, motorCyclePoint1, edgePoint0, edgePoint1, out float motorcycleTime, out float edgeTime))
				{
					// we hit the actual edge
					if (0.0f <= edgeTime && edgeTime <= 1.0f)
					{
						// don't hit adjacent edges that we originate from - so we make sure that we have at least moved a LITTLE forward
						if (motorcycleTime > 0.001f)
						{
							motorcycleTime += motorcycle.startTime;
							if (motorcycleTime < motorcycle.crashTime)
							{
								motorcycle.UpdateCrash(MotorcycleCrashType.Wall, i, motorcycleTime);
							}
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
					if (crash.time < motorcycle.crashTime)
					{
						motorcycle.UpdateCrash(MotorcycleCrashType.Motorcycle, crash.motorcycleTraceIndex, crash.time);
					}

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
						nUnprocessedMotorcycles += ProcessCrashedMotorcycles(realCrashes, newMotorcycles, nSimultaneousCrashes, i);
						nSimultaneousCrashes = 1;
					}
				}

				// we process the last batch of motorcycles
				nUnprocessedMotorcycles += ProcessCrashedMotorcycles(realCrashes, newMotorcycles, nSimultaneousCrashes, realCrashes.Count - 1);

				// append the new motorcycles to the old ones
				Motorcycle[] allMotorcycles = new Motorcycle[nMotorcycles + nUnprocessedMotorcycles];
				Array.Copy(motorcycles, allMotorcycles, nMotorcycles);
				for (int i = 0; i < newMotorcycles.Count; ++i)
				{
					allMotorcycles[nMotorcycles + i] = newMotorcycles[i];

					// make sure we also calculate wall collisions for these new motorcycles!
					ref Motorcycle motorcycle = ref allMotorcycles[nMotorcycles + i];
					CollideMotorcycleWithPolygon(ref motorcycle);
				}

				motorcycles = allMotorcycles;
				nMotorcycles = motorcycles.Length;
			}

			// finally, we go over all motorcycles that crashed against a WALL, and create a new vertex at the wall - this might split up an edge!
			// we already processed motorcycles that crashed against each other or a trace of another motorcycle
			for (int i = 0; i < motorcycles.Length; ++i)
			{
				ref var motorcycle = ref motorcycles[i];
				if (motorcycle.crashType == MotorcycleCrashType.Wall)
				{
					SplitEdgeAtCrashSite(motorcycle);
				}
			}
		}

		private int ProcessCrashedMotorcycles(List<MotorcycleCrash> realCrashes, List<Motorcycle> newMotorcycles, int nSimultaneousCrashes, int lastMotorcycleIndex)
		{
			// we add new edges for the trace of each motorcycle
			int newVertexIndex = AddNewEdges(realCrashes, lastMotorcycleIndex - nSimultaneousCrashes + 1, nSimultaneousCrashes);

			int nUnprocessedMotorcycles = 0;
			if (nSimultaneousCrashes > 1)
			{
				nUnprocessedMotorcycles = AddNewMotorcycles(realCrashes, lastMotorcycleIndex - nSimultaneousCrashes + 1, nSimultaneousCrashes, newVertexIndex, newMotorcycles);
			}

			return nUnprocessedMotorcycles;
		}

		private int AddNewEdges(List<MotorcycleCrash> crashes, int startIndex, int nSimultaneousCrashes)
		{
			// firstly, we add a new vertex at the crash point - might be shared by multiple edges!
			// velocity will be calculated later
			MotorcycleCrash crash = crashes[startIndex];
			WavefrontVertex crashVertex = new WavefrontVertex(crash.crashPos, float2.zero, WavefrontVertexType.SteinerResting);
			vertices.Add(crashVertex);
			int crashVertexIndex = vertices.Count - 1;

			// we add a new edge for each crashed motorcycle
			for (int i = startIndex; i < startIndex + nSimultaneousCrashes; ++i)
			{
				ref var motorcycle = ref motorcycles[crashes[i].motorcycleIndex];
				WavefrontEdge newEdge = new WavefrontEdge(motorcycle.vertexIndex, crashVertexIndex);
				edges.Add(newEdge);
			}

			return crashVertexIndex;
		}

		private void SplitEdgeAtCrashSite(in Motorcycle motorcycle)
		{
			if (motorcycle.crashType != MotorcycleCrashType.Wall) throw new ArgumentException($"Motorcycle {motorcycle} did not crash against a wall.");

			WavefrontEdge crashEdge = edges[motorcycle.crashTargetIndex];
			float2 crashPos = motorcycle.getCrashPos();
			int startVertex = motorcycle.vertexIndex;

			// get the start and endpoints of the crash edge
			float2 p0 = vertices[crashEdge.prevVertexIndex].pos;
			float2 p1 = vertices[crashEdge.nextVertexIndex].pos;

			// if we crashed right into an endpoint of the edge, we don't need to create a new vertex
			if (math.distance(p0, crashPos) < Geometry.EPS)
			{
				edges.Add(new WavefrontEdge(startVertex, crashEdge.prevVertexIndex));
			}
			else if (math.distance(p1, crashPos) < Geometry.EPS)
			{
				edges.Add(new WavefrontEdge(startVertex, crashEdge.nextVertexIndex));
			}

			// this is the complicated part - we need to split up the edge!
			else
			{
				WavefrontVertex newCrashVertex = new WavefrontVertex(crashPos, float2.zero, WavefrontVertexType.SteinerMoving);
				vertices.Add(newCrashVertex);
				int newVertexIndex = vertices.Count - 1;

				// first, add the motorcycle trace edge
				edges.Add(new WavefrontEdge(startVertex, newVertexIndex));

				// now, we repurpose the original edge (so we don't need to fuck with updating indices) to be one part of the new edge
				edges[motorcycle.crashTargetIndex] = new WavefrontEdge(crashEdge.prevVertexIndex, newVertexIndex);
				edges.Add(new WavefrontEdge(newVertexIndex, crashEdge.nextVertexIndex));
			}
		}

		private List<Vector2> incomingMotorcycleDirs = new List<Vector2>();
		private int AddNewMotorcycles(List<MotorcycleCrash> crashes, int startIndex, int nSimultaneousCrashes, int startVertexIndex, List<Motorcycle> newMotorcycles)
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
			for (int i = 0; i < incomingMotorcycleDirs.Count; ++i)
			{
				float2 dir0 = incomingMotorcycleDirs[i];
				float2 dir1 = incomingMotorcycleDirs[(i + 1) % incomingMotorcycleDirs.Count];

				float2 prev = p - dir1;
				float2 curr = p;
				float2 next = p - dir0;

				// we flip them because we need to find the relfex 
				if (Geometry.IsRelfexVertex(prev, curr, next))
				{
					float2 newVelocity = Wavefront.CalculateVelocity(prev, curr, next);
					var newMotorcycle = new Motorcycle(startVertexIndex, crashTime, p, newVelocity);
					newMotorcycles.Add(newMotorcycle);
					++nNewMotorcycles;
				}
			}

			return nNewMotorcycles;

		}
	}
}