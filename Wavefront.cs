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
	public class Wavefront
	{
		public readonly WavefrontVertex[] vertices;


		public Wavefront(float2[] contour)
		{

			// convert all vertices in the contour to wavefront vertices
			vertices = new WavefrontVertex[contour.Length];
			for (int i = 0; i < contour.Length; ++i)
			{
				float2 prev = contour[(i - 1 + contour.Length) % contour.Length];
				float2 curr = contour[i];
				float2 next = contour[(i + 1) % contour.Length];

				/*float2 prevToCurr = curr - prev;
				float2 nextToCurr = next - curr;

				// rotate them 90°
				float2 moveDir1 = math.normalize(Geometry.RotateMinus90Degrees(prevToCurr));
				float2 moveDir2 = math.normalize(Geometry.RotateMinus90Degrees(nextToCurr));

				float2 velocity = moveDir1 + moveDir2;
				float speed = math.length(velocity);*/

				float2 velocity = CalculateVelocity(prev, curr, next);

				bool isReflex = Geometry.IsRelfexVertex(prev, curr, next);
				vertices[i] = new WavefrontVertex(curr, velocity, isReflex ? WavefrontVertexType.Reflex : WavefrontVertexType.Convex);
			}
		}

		public static float2 CalculateVelocity(in float2 prev, in float2 curr, in float2 next)
		{
			float2 prevToCurr = curr - prev;
			float2 nextToCurr = next - curr;

			// rotate them 90°
			float2 moveDir1 = math.normalize(Geometry.RotateMinus90Degrees(prevToCurr));
			float2 moveDir2 = math.normalize(Geometry.RotateMinus90Degrees(nextToCurr));

			float2 dir = math.normalize(moveDir1 + moveDir2);
			float speed = 2.0f / (math.dot(moveDir1, dir) + math.dot(moveDir2, dir));
			float2 velocity = dir * speed;

			return velocity;
		}
	}
}