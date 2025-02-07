using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Profiling;
using UnityEngine.Rendering;
using UnityEngine.UI;

namespace Briganti.StraightSkeletons
{
	public class Geometry
	{
		public const float EPS = 0.0001f;

		public static bool GetLineIntersection(float2 p0, float2 p1, float2 q0, float2 q1, out float t0, out float t1)
		{
			// http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

			// calculate the direction vectors
			float2 r = p1 - p0;
			float2 s = q1 - q0;

			// calculate the 2D cross product between the dir vectors
			float dirCross = Cross2D(r, s);

			// no intersection - they are colinear or parallel
			if (math.abs(dirCross) < math.EPSILON)
			{
				t0 = float.PositiveInfinity;
				t1 = float.PositiveInfinity;
				return false;
			}

			// calculate the u-value
			t0 = Cross2D(q0 - p0, s) / dirCross;
			t1 = Cross2D(q0 - p0, r) / dirCross;
			return true;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public static float Cross2D(float2 v, float2 w)
		{
			return v.x * w.y - v.y * w.x;
		}

		// returns t0 for the line p0 -> p1
		public static bool GetLineSegmentIntersection(float2 p0, float2 p1, float2 q0, float2 q1, ref float t0, ref float t1)
		{
			GetLineIntersection(p0, p1, q0, q1, out t0, out t1);
			if (t0 < 0.0f || t0 > 1.0f || t1 < 0.0f || t1 > 1.0f) return false;
			return true;
		}

		public static bool IsCounterClockwise2D(IList<float2> points)
		{
			// https://www.element84.com/blog/determining-the-winding-of-a-polygon-given-as-a-set-of-ordered-points
			// https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
			float sum = 0.0f;
			for (int i = 0; i < points.Count; ++i)
			{
				sum += (points[(i + 1) % points.Count].x - points[i].x) * (points[(i + 1) % points.Count].y + points[i].y);
			}

			return sum < 0;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public static float Det(float2 prev, float2 curr, float2 next)
		{
			float2 a = prev;
			float2 b = curr;
			float2 c = next;
			return (b.x - a.x) * (c.y - b.y) - (c.x - b.x) * (b.y - a.y);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public static bool IsRelfexVertex(float2 prev, float2 curr, float2 next)
		{
			return Det(prev, curr, next) > 0;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public static float2 Rotate90Degrees(float2 v)
		{
			return new float2(-v.y, v.x);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public static float2 RotateMinus90Degrees(float2 v)
		{
			return new float2(v.y, -v.x);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public static float2 Rotate90Degrees(float2 v, bool positiveDirection)
		{
			if (positiveDirection)
			{
				return new float2(-v.y, v.x);
			}
			else
			{
				return new float2(v.y, -v.x);
			}
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public static float GetAngle(float2 v)
		{
			return math.atan2(v.y, v.x);
		}

		// project a point on a line
		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public static float2 ProjectPointOnLine(float2 p1, float2 p2, float2 p, out float t)
		{
			// http://paulbourke.net/geometry/pointline/

			float2 p3 = p;
			float u = ((p3.x - p1.x) * (p2.x - p1.x) + (p3.y - p1.y) * (p2.y - p1.y)) / math.lengthsq(p2 - p1);
			float2 proj = new(p1.x + u * (p2.x - p1.x), p1.y + u * (p2.y - p1.y));
			t = u;
			return proj;
		}

		public static int2 RoundToInt(float2 v)
		{
			return new int2((int)math.round(v.x), (int)math.round(v.y));
		}
	}
}