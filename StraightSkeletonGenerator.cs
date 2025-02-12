
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
		private VertexGraph graph;

		private FastPriorityQueue<Edge> priorityQueue;

		private float time = 0.0f;


		public StraightSkeletonGenerator(PolygonWithHoles polygonWithHoles)
		{

			// first, convert into a polygon wavefront
			PolygonWavefront polygonWavefront = new PolygonWavefront(polygonWithHoles);

			// then, calculate the motorcycles
			MotorcycleGraphGenerator motorcycleGraphGenerator = new MotorcycleGraphGenerator(polygonWavefront);
			graph = motorcycleGraphGenerator.CalculateWavefrontWithMotorcycles();

			// set up the priority queue - we never have more edges than at the start
			priorityQueue = new FastPriorityQueue<Edge>(graph.edges.Count);

			// we update the collapse time & add them to the priority queue
			for (int i = 0; i < graph.edges.Count; ++i)
			{
				UpdateCollapseTime(i);
				priorityQueue.Enqueue(graph.edges[i], graph.edges[i].collapseTime);
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
			Edge edge = graph.edges[edgeIndex];

			int prevVertexIndex = edge.prevVertexIndex;
			int nextVertexIndeg = edge.nextVertexIndex;
			ref var prevVertex = ref graph.vertices[edge.prevVertexIndex];
			ref var nextVertex = ref graph.vertices[edge.nextVertexIndex];

			// if both are NOT wavefront vertices, they are not moving
			if (!prevVertex.IsWavefrontVertex() && !nextVertex.IsWavefrontVertex())
			{
				edge.UpdateCollapse(float.MaxValue, float2.zero);
				return;
			}

			// we know both velocities - this is a simple edge event and we can easily calculate what's going to happen
			if (prevVertex.HasVelocity() && nextVertex.HasVelocity())
			{
				GetCollapseBetweenTwoVelocities(prevVertex, nextVertex, out float collapseTime, out float2 collapsePos);
				edge.UpdateCollapse(collapseTime, collapsePos);
				return;
			}

			// one of the vertices does NOT have a velocity - we can still easily calculate the collapse point based on the intersection of the dir and the velocity
			// this is still a case where we are on the wavefront edge, since at least one has velocity
			// this is an edge between a convex/reflex vertex and a moving steiner vertex
			if (nextVertex.HasVelocity())
			{

				// the other side is a moving steiner vertex - in that case, we 

				/*GetCollapseBetweenVelocityAndDir(prevVertexDir, nextVertex, out float collapseTime, out float2 collapsePos);
				edge.UpdateCollapse(collapseTime, collapsePos);*/
				return;
			}

			// if any of the vertices is resting, the edge isn't moving
			edge.UpdateCollapse(float.MaxValue, float2.zero);
			if (!prevVertex.IsMoving() && !nextVertex.IsMoving()) return;

			// this does NOT work since the velocity of the vertex is different for each adjacent edge???

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

		private void CalculateDirandSpeed(ref Edge edge, int vertexIndex, out float2 dir, out float speed) {
			ref Vertex vertex = ref graph.vertices[vertexIndex];

			// if this vertex has a velocity, it is an original polygon vertex - we just return the velocity!
			if (vertex.HasVelocity())
			{
				dir = math.normalize(vertex.velocity);
				speed = math.length(vertex.velocity);
			}

			// otherwise, we are dealing with a moving steiner vertex, and we move along the edge of the motorcycle graph at an unknown velocity
			if (vertex.type != WavefrontVertexType.SteinerMoving) throw new ArgumentException($"Cannot calculate the dir/speed of vertex {vertex} because it is not an original wavefront vertex or a moving steiner vertex.");
			dir = graph.GetMotorcycleDir(edge, vertexIndex);
			speed = 0.0f;
		}

		private void GetCollapseBetweenTwoVelocities(in Vertex prevVertex, in Vertex nextVertex, out float collapseTime, out float2 collapsePos)
		{
			collapseTime = float.MaxValue;
			collapsePos = float2.zero;

			if (Geometry.GetLineIntersection(prevVertex.pos, prevVertex.pos + prevVertex.velocity, nextVertex.pos, nextVertex.pos + nextVertex.velocity, out float t0, out float t1))
			{
				// our intersection is behind us - we don't collapse
				if (t0 <= 0.0f || t1 <= 0.0f) return;

				// we calculate the collapse time
				collapseTime = time + t0;
				collapsePos = prevVertex.pos + prevVertex.velocity * t0;
			}
		}

		/*private void GetCollapseBetweenVelocityAndDir(Edge edge, in Vertex steinerVertex, in Vertex normalVertex, out float collapseTime, out float2 collapsePos)
		{
			collapseTime = float.MaxValue;
			collapsePos = float2.zero;

			if (Geometry.GetLineIntersection(prevVertex.pos, prevVertex.pos + prevVertex.velocity, nextVertex.pos, nextVertex.pos + nextVertex.velocity, out float t0, out float t1))
			{
				// our intersection is behind us - we don't collapse
				if (t0 <= 0.0f || t1 <= 0.0f) return;

				// we calculate the collapse time
				collapseTime = time + t0;
				collapsePos = prevVertex.pos + prevVertex.velocity * t0;
			}
		}

		private float2 GetVertexVelocity(Edge edge, int vertexIndex)
		{
			ref Vertex vertex = ref graph.vertices[vertexIndex];

			switch (vertex.type)
			{

				// if this is an original polygon vertex, we just have a velocity - it's easy, we already have a velocity defined
				case WavefrontVertexType.Convex:
				case WavefrontVertexType.Reflex:
					return vertex.velocity;

				case WavefrontVertexType.ConvexAndSteiner:
					// TODO
					throw new NotImplementedException();

				// moving steiner vertices can have a different velocity depending on the edge that is adjacent to the edge we're trying to calculate
				case WavefrontVertexType.SteinerMoving:

				// these are just not moving right now
				case WavefrontVertexType.SteinerResting:
				case WavefrontVertexType.SteinerMulti:
					return float2.zero;

				default:
					throw new NotSupportedException($"Vertex type {vertex.type} is not supported.");

			}
		}

		private float2 GetMovingSteinerVelocity(Edge edge, int vertexIndex)
		{
			ref Vertex vertex = ref graph.vertices[vertexIndex];
			if (vertex.type != WavefrontVertexType.SteinerMoving) throw new ArgumentException($"Vertex {vertex} is not a moving steiner vertex!");

			// get the direction in which the steiner vertex moves - this is the adjacent steiner edge
			Edge adjacentEdge = edge.GetAdjacentEdge(vertexIndex);
			if (graph.IsWavefrontEdge(adjacentEdge)) throw new ArgumentException($"Adjacent edge {adjacentEdge} should NOT be a wavefront edge but a steiner edge!");
			float2 dir = graph.GetDirFrom(adjacentEdge, vertexIndex);

			// if the edge we're calculating the velocity for is a wavefront edge, we know that the wavefront is moving perpendicular to the edge
			// this makes it easy to calculate our own velocity
			if (graph.IsWavefrontEdge(edge))
			{
				float2 wavefrontDir = graph.GetWavefrontDir(edge);

				// use cas from soscastoa to calculate the actual velocity of the moving steiner vertex
				float angle = Geometry.Angle(wavefrontDir, dir);
				float2 velocity = math.cos(angle) * 1.0f; // the wavefront moves at speed 1
				return velocity;
			}

			// ... however, IF the edge we're calculating the velocity for is a steiner edge,
			// we need to calculate the distance it takes to move to the next resting/multi steiner vertex




			// now we calculate the actual position of the wavefront emanating from this vertex

			return

			edge.prevClockwiseEdge;
		}

		private float GetMovingSteinerSpeedAlongSteinerEdge(Edge edge, int vertexIndex)
		{
			if (graph.IsWavefrontEdge(edge)) throw new ArgumentException($"Edge {edge} is not an edge between two steiner vertices.");

			ref Vertex vertex = ref graph.vertices[vertexIndex];

			// this vertex is not moving yet
			if (!vertex.IsMoving()) return 0;

			// if it IS moving, it MUST be moving along this edge
			float2 dir = graph.GetDirFrom(edge, vertexIndex);

			// but we also need to know the speed at which the wavefront is moving at this vertex
		}*/

		/*private float2 GetMoveDirAtVertex(Edge edge, int vertexIndex)
		{
			bool isWavefrontEdge = IsWavefrontEdge(edge);

			ref Vertex vertex = ref vertices[vertexIndex];
			if (vertex.type != WavefrontVertexType.SteinerMoving)
			{
				return math.normalize(vertex.velocity);
			}
			else
			{

			}
		}*/

		private void ProcessEvent(Edge edge)
		{
			/*

			// we update the current time to the collapse time of this edge
			time = edge.collapseTime;

			var prevVertex = vertices[edge.prevVertexIndex];
			var nextVertex = vertices[edge.nextVertexIndex];

			// 1. edge event - two simple convex vertices of the original wavefront meat each other and the edge between them disappears
			if (prevVertex.type == WavefrontVertexType.Convex && nextVertex.type == WavefrontVertexType.Convex)
			{
				ProcessClassicEdgeEvent(edge);
			}*/
		}

	}
}