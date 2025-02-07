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
	public enum WavefrontVertexType
	{
		Convex,
		Reflex,
		ConvexMulti,
		SteinerMoving,
		SteinerResting,
		SteinerMulti
	}

	public static class WavefrontVertexTypeExtensions
	{
		public static bool IsMoving(this WavefrontVertexType type)
		{
			return type != WavefrontVertexType.SteinerResting && type != WavefrontVertexType.SteinerMulti;
		}
	}
}