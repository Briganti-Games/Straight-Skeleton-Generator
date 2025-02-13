using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Profiling;
using UnityEngine.Rendering;
using UnityEngine.UI;

namespace Briganti.StraightSkeletonGeneration
{
	public class PolygonWithHoles
	{
		public float2[] outerContour;
		public List<float2[]> innerContours = new List<float2[]>();
	}
}