using System.Collections.Generic;
using Unity.Mathematics;

namespace Briganti.StraightSkeletonGeneration
{
	public class PolygonWithHoles
	{
		public float2[] outerContourCounterClockwise;
		public List<float2[]> innerContoursClockwise = new List<float2[]>();

		// set everything manually
		public PolygonWithHoles()
		{
		}

		// only an outer contour
		public PolygonWithHoles(float2[] outerContourCounterClockwise)
		{
			this.outerContourCounterClockwise = outerContourCounterClockwise;
		}

		// outer + inner contours
		public PolygonWithHoles(float2[] outerContour, List<float2[]> innerContoursClockwise)
		{
			this.outerContourCounterClockwise = outerContour;
			this.innerContoursClockwise = innerContoursClockwise;
		}
	}
}