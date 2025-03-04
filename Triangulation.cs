using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using UnityEngine.Profiling;

namespace Briganti.StraightSkeletonGeneration
{
	public static class Triangulation
	{

		public static int[] Triangulate(Polygon polygon)
		{
			if (polygon.IsHorizontallyMonotone())
			{
				return TriangulateHorizontallyMonotone(polygon);
			}
			else
			{
				throw new NotImplementedException("Cannot triangulate non-monotone polygon: must implement Ear Clipping or some other algorithm for these cases...");
			}

		}

		private static int[] TriangulateHorizontallyMonotone(Polygon polygon)
		{
			Profiler.BeginSample("Triangulation.TriangulateHorizontallyMonotone");

			// https://www.cs.umd.edu/class/spring2020/cmsc754/Lects/lect05-triangulate.pdf
			// https://sites.cs.ucsb.edu/~suri/cs235/Triangulation.pdf

			// in this implementation, we assume Polygon is counter-clockwise, which is enforced by the Polygon class so should be safe!

			// first, we find the leftmost vertex in the polygon
			int leftMostIndex = 0;
			float minX = float.MaxValue;
			var vertices = polygon.vertices;
			int n = vertices.Count;
			if (n < 3) throw new ArgumentException($"Cannot triangulate a polygon with less than 3 vertices!");
			for (int i = 0; i < n; ++i)
			{
				float x = vertices[i].x;
				if (x < minX)
				{
					minX = x;
					leftMostIndex = i;
				}
			}

			// now we have two chains, (bottom and top), and we traverse them step by step and create triangles
			List<int> triangles = new List<int>(vertices.Count*3);
			int nextBottomIndex = (leftMostIndex + 1) % n;
			int nextTopIndex = (leftMostIndex - 1 + n) % n;

			// a "reflex" vertex is a vertex which points "inwards" towards the polygon
			// this can be determined by seeing if the cross product of the point with its surrounding points is counter-clockwise
			List<int> reflexChain = new List<int>(vertices.Count);
			reflexChain.Add(leftMostIndex);
			bool reflexIsTopChain = true;

			bool IsCounterClockwise(int i, int j, int k)
			{
				Vector2 a = vertices[i];
				Vector2 b = vertices[j];
				Vector2 c = vertices[k];

				return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y) > 0;
			}

			for (int i = 1; i < n; ++i)
			{
				bool topChain;
				int nextIndex;

				// these offsets are used to invert the triangle order when doing a bottom chain connection,
				// because we need to ensure all final triangles are counter-clockwise!
				int offset0 = 0, offset1 = 0;

				if (vertices[nextBottomIndex].x <= vertices[nextTopIndex].x)
				{
					topChain = false;
					nextIndex = nextBottomIndex;
					nextBottomIndex = (nextBottomIndex + 1) % n;
					offset0 = 1;
					offset1 = -1;
				}
				else
				{
					topChain = true;
					nextIndex = nextTopIndex;
					nextTopIndex = (nextTopIndex - 1 + n) % n;
				}

				// if this is the second one, we just add it and are done
				// we need 2 vertices to get the party started
				if (i == 1)
				{
					reflexIsTopChain = topChain;
					reflexChain.Add(nextIndex);
				}

				// otherwise, we see if we need to connect
				else
				{
					if (topChain == reflexIsTopChain)
					{

						// we need to connect vertices until we hit a clockwise vertex
						while (reflexChain.Count > 1 && IsCounterClockwise(reflexChain[reflexChain.Count - 2 + offset0], nextIndex, reflexChain[reflexChain.Count - 1 + offset1]))
						{
							triangles.Add(reflexChain[reflexChain.Count - 2 + offset0]);
							triangles.Add(nextIndex);
							triangles.Add(reflexChain[reflexChain.Count - 1 + offset1]);
							reflexChain.RemoveAt(reflexChain.Count - 1);
						}
						reflexChain.Add(nextIndex);
					}

					// we are on the other side of the chain - we connect ALL vertices in the chain and clear the chain
					else
					{

						for (int j = 0; j < reflexChain.Count - 1; ++j)
						{
							triangles.Add(reflexChain[j + offset0]);
							triangles.Add(reflexChain[j + 1 + offset1]);
							triangles.Add(nextIndex);
						}

						int lastReflexChainIndex = reflexChain[reflexChain.Count - 1];
						reflexChain.Clear();
						reflexChain.Add(lastReflexChainIndex);
						reflexChain.Add(nextIndex);

						reflexIsTopChain = topChain;
					}
				}
			}

			Profiler.EndSample();
			return triangles.ToArray();
		}
	}
}