
using Briganti.StraightSkeletons.Priority_Queue;
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
		private Wavefront wavefront;
		private List<WavefrontVertex> vertices => wavefront.vertices;
		private List<WavefrontEdge> edges => wavefront.edges;

		private FastPriorityQueue<WavefrontEdge> priorityQueue;

		private float time = 0.0f;


		public StraightSkeletonGenerator(PolygonWithHoles polygonWithHoles)
		{

			// first, convert into a polygon wavefront
			PolygonWavefront polygonWavefront = new PolygonWavefront(polygonWithHoles);

			// then, calculate the motorcycles
			MotorcycleGraphGenerator motorcycleGraphGenerator = new MotorcycleGraphGenerator(polygonWavefront);
			wavefront = motorcycleGraphGenerator.CalculateWavefrontWithMotorcycles();

			// set up the priority queue - we never have more edges than at the start
			priorityQueue = new FastPriorityQueue<WavefrontEdge>(edges.Count);

			// we update the collapse time & add them to the priority queue
			for (int i = 0; i < edges.Count; ++i)
			{
				UpdateCollapseTime(i);
				priorityQueue.Enqueue(edges[i], edges[i].collapseTime);
			}

			// now we keep going until the queue is empty
			while (priorityQueue.Count > 0)
			{
				var edge = priorityQueue.Dequeue();
				ProcessEvent(edge);
			}
		}

		public void UpdateCollapseTime(int edgeIndex)
		{
			var edge = edges[edgeIndex];
			var prevVertex = vertices[edge.prevVertexIndex];
			var nextVertex = vertices[edge.nextVertexIndex];

			// if any of the vertices is resting, the edge isn't moving
			edge.UpdateCollapse(float.MinValue, float2.zero);
			if (!prevVertex.type.IsMoving() && !nextVertex.type.IsMoving()) return;

			// calculate the intersection point of the two vertex movements
			if (math.length(prevVertex.velocity) < Geometry.EPS || math.length(nextVertex.velocity) < Geometry.EPS) throw new ArgumentException($"Edge {edge} does not have two moving vertices: {prevVertex} and {nextVertex}");

			// we calculate the intersection point of the two velocities
			if (Geometry.GetLineIntersection(prevVertex.pos, prevVertex.pos + prevVertex.velocity, nextVertex.pos, nextVertex.pos + nextVertex.velocity, out float t0, out float t1))
			{
				if (t0 <= 0.0f || t1 <= 0.0f) return;

				// we calculate the collapse time
				float collapseTime = time + t0;
				float2 collapsePos = prevVertex.pos + prevVertex.velocity * t0;
				edge.UpdateCollapse(collapseTime, collapsePos);
			}
		}

		private void ProcessEvent(WavefrontEdge edge) {


			// we update the current time to the collapse time of this edge
			time = edge.collapseTime;

			var prevVertex = vertices[edge.prevVertexIndex];
			var nextVertex = vertices[edge.nextVertexIndex];

			// 1. edge event
			if (prevVertex.type == WavefrontVertexType.Convex && nextVertex.type == WavefrontVertexType.Convex)
			{
				ProcessEdgeEvent();
			}
		}

		private void ProcessEdgeEvent() {

		}
	}
}