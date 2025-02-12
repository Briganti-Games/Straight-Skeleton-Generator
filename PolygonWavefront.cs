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
	public class PolygonWavefront
	{
		public readonly PolygonWithHoles polygonWithHoles;

		public Vertex[] vertices;


		public List<VertexGraph> wavefronts = new List<VertexGraph>();


		public PolygonWavefront(PolygonWithHoles polygonWithHoles)
		{
			this.polygonWithHoles = polygonWithHoles;

			// make sure the outer contour is clockwise and the inner oners counter clockwise
			if (Geometry.IsCounterClockwise2D(polygonWithHoles.outerContour)) throw new ArgumentException($"The outer contour is not defined clockwise!");
			wavefronts.Add(new VertexGraph(polygonWithHoles.outerContour));

			for (int i = 0; i < polygonWithHoles.innerContours.Count; ++i)
			{
				var innerContour = polygonWithHoles.innerContours[i];
				if (!Geometry.IsCounterClockwise2D(innerContour)) throw new ArgumentException($"Inner contour {i} is not defined counter-clockwise.");
				wavefronts.Add(new VertexGraph(polygonWithHoles.innerContours[i]));
			}
		}
	}
}