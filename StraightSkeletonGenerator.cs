
using System;
using System.Collections.Generic;
using System.Linq;
using Unity.Collections;
using Unity.Mathematics;
using UnityEngine;

namespace Briganti.StraightSkeletons
{
	public class StraightSkeletonGenerator
	{
		public StraightSkeletonGenerator(PolygonWithHoles polygonWithHoles)
		{

			// first, convert into a polygon wavefront
			PolygonWavefront polygonWavefront = new PolygonWavefront(polygonWithHoles);

			// then, calculate the motorcycles
			MotorcycleGraphGenerator motorcycleGraphGenerator = new MotorcycleGraphGenerator(polygonWavefront);
			Motorcycle[] motorcycles = motorcycleGraphGenerator.CalculateGraph();

			// we now convert the polygon wavefront into a list of wavefronts (polygons) that are all convex
			var wavefronts = polygonWavefront.MergeMotorcycles(motorcycles);
		}
	}
}