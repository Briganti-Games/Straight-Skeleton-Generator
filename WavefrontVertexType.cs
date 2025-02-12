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
		ConvexAndSteiner,
		SteinerMoving,
		SteinerResting,
		SteinerMulti
	}

	public static class WavefrontVertexTypeExtensions
	{
	}
}