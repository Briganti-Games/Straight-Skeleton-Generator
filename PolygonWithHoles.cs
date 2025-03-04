using System.Collections.Generic;
using Unity.Mathematics;

namespace Briganti.StraightSkeletonGeneration
{
	public class PolygonWithHoles
	{
		public float2[] outerContour;
		public List<float2[]> innerContours = new List<float2[]>();
	}
}