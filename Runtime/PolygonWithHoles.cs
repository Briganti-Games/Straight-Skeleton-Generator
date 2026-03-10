using System.Collections.Generic;
using Unity.Mathematics;

namespace Briganti.StraightSkeletonGeneration
{
	public class PolygonWithHoles
	{
		public float2[] outerContourClockwise;
		public List<float2[]> innerContoursCounterClockwise = new List<float2[]>();

		// set everything manually
		public PolygonWithHoles()
		{
		}

		// only an outer contour
		public PolygonWithHoles(float2[] outerContourClockwise)
		{
			this.outerContourClockwise = outerContourClockwise;
		}

		// outer + inner contours
		public PolygonWithHoles(float2[] outerContourClockwise, List<float2[]> innerContoursCounterClockwise)
		{
			this.outerContourClockwise = outerContourClockwise;
			this.innerContoursCounterClockwise = innerContoursCounterClockwise;
		}
	}
}