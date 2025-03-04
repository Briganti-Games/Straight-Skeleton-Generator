using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using Unity.Mathematics;
using UnityEngine;

namespace Briganti.StraightSkeletonGeneration
{
	public class StraightSkeletonMeshGenerator
	{
		private readonly StraightSkeleton straightSkeleton;

		private int maxEdgesStartingFromVertex = 1;

		public class Slab
		{
			public List<int> indices = new List<int>();
		}


		public StraightSkeletonMeshGenerator(StraightSkeleton straightSkeleton)
		{
			this.straightSkeleton = straightSkeleton;

			// calculate the max vertex degree - this is used to create a convenient mapping of vertex to edges
			int[] vertexDegrees = new int[straightSkeleton.nVertices];
			for (int i = 0; i < straightSkeleton.nEdges; ++i)
			{
				ref var edge = ref straightSkeleton.edges[i];
				vertexDegrees[edge.prevVertexIndex]++;
			}

			maxEdgesStartingFromVertex = vertexDegrees.Max();
		}

		public Mesh GenerateRoofMesh(float maxTimeHeight = 1, float texCoordScale = 1)
		{
			List<Vector3> meshVertices = new List<Vector3>();
			List<Vector2> meshUVs = new List<Vector2>();
			List<int> meshTriangles = new List<int>();
			List<Vector3> meshNormals = new List<Vector3>();

			void AddTriangle(int v1, int v2, int v3)
			{
				meshTriangles.Add(v1);
				meshTriangles.Add(v2);
				meshTriangles.Add(v3);
			}

			List<Slab> slabs = GenerateStraightSkeletonSlabs();

			// first, we need to triangulate the different slabs
			for (int slabIndex = 0; slabIndex < slabs.Count; ++slabIndex)
			{
				int startMeshIndex = meshVertices.Count;
				Slab slab = slabs[slabIndex];
				if (slab.indices.Count < 3) throw new ArgumentException($"Slab {string.Join(" -> ", slab.indices)} does not the minimum of 3 vertices to form a triangle mesh.");

				float2 origin = straightSkeleton.vertices[slab.indices[0]];
				float2 dir = math.normalize(origin - straightSkeleton.vertices[slab.indices[slab.indices.Count - 1]]);
				float2 perp = new float2(dir.y, -dir.x);

				List<Vector2> slabVertices = new List<Vector2>(slab.indices.Count);
				for (int i = 0; i < slab.indices.Count; ++i)
				{
					int vertexIndex = slab.indices[i];

					float2 vertex = straightSkeleton.vertices[vertexIndex];
					float time = straightSkeleton.vertexTimes[vertexIndex];

					Vector2 slabVertex = new Vector2(vertex.x, vertex.y);
					slabVertices.Add(slabVertex);

					float y = time / straightSkeleton.maxTime * maxTimeHeight;
					Vector3 meshVertex = new Vector3(vertex.x, y, vertex.y);
					meshVertices.Add(meshVertex);

					float2 relative = vertex - origin;
					meshUVs.Add(new Vector2(math.dot(relative, dir), math.dot(relative, perp) * texCoordScale));
				}

				// calculate the normal for this slab, which is uniform along all its vertices
				int v1 = startMeshIndex;
				int v2 = startMeshIndex + 1;
				int v3 = startMeshIndex + 2;
				Vector3 normal = Vector3.Cross(meshVertices[v2] - meshVertices[v1], meshVertices[v3] - meshVertices[v1]);
				for (int i = 0; i < slab.indices.Count; ++i)
				{
					meshNormals.Add(normal);
				}

				try
				{
					Polygon polygon = new Polygon(slabVertices);
					int[] triangleIndices = Triangulation.Triangulate(polygon);
					for (int i = 0; i < triangleIndices.Length; i += 3)
					{
						AddTriangle(startMeshIndex + triangleIndices[i], startMeshIndex + triangleIndices[i + 1], startMeshIndex + triangleIndices[i + 2]);
					}
				}
				catch (NotImplementedException e)
				{
					// we failed to triangulate - we print an error but we proceed
					Debug.Log($"Failed to triangulate slab {string.Join(" -> ", slabVertices)} with indices {string.Join(" -> ", slab.indices)}: {e.Message}");
				}
			}

			Mesh mesh = new Mesh();
			mesh.SetVertices(meshVertices);
			mesh.SetNormals(meshNormals);
			mesh.SetUVs(0, meshUVs);
			mesh.SetTriangles(meshTriangles, 0);
			return mesh;
		}

		public List<Slab> GenerateStraightSkeletonSlabs()
		{
			List<Slab> polygons = new List<Slab>();

			int[] vertexEdgeMapping = new int[maxEdgesStartingFromVertex * straightSkeleton.nVertices];

			// reset
			for (int i = 0; i < vertexEdgeMapping.Length; ++i)
			{
				vertexEdgeMapping[i] = -1;
			}

			// map every edge to the right vertices
			bool[] edgesVisited = new bool[straightSkeleton.nEdges];
			for (int edgeIndex = 0; edgeIndex < straightSkeleton.nEdges; ++edgeIndex)
			{
				ref var edge = ref straightSkeleton.edges[edgeIndex];
				int prevVertexIndex = edge.prevVertexIndex;

				int idx = prevVertexIndex * maxEdgesStartingFromVertex;
				while (vertexEdgeMapping[idx] != -1) ++idx;
				if (idx >= (prevVertexIndex + 1) * maxEdgesStartingFromVertex) throw new InvalidOperationException($"There is a vertex {prevVertexIndex} with more incoming edges than {maxEdgesStartingFromVertex}!");
				vertexEdgeMapping[idx] = edgeIndex;
			}

			int GetNextEdge(int edgeIndex)
			{
				ref var edge = ref straightSkeleton.edges[edgeIndex];
				float2 edgeDir = math.normalize(straightSkeleton.vertices[edge.nextVertexIndex] - straightSkeleton.vertices[edge.prevVertexIndex]);

				int nextVertexIndex = edge.nextVertexIndex;
				int idx = nextVertexIndex * maxEdgesStartingFromVertex;

				// go over all edges connected to this vertex, and find the next one clockwise from this one
				int nextEdgeIndex = -1;
				float bestAngle = float.MinValue;
				while (idx < (nextVertexIndex + 1) * maxEdgesStartingFromVertex && vertexEdgeMapping[idx] != -1)
				{
					int nextEdgeCandidateIndex = vertexEdgeMapping[idx];
					++idx;

					if (nextEdgeCandidateIndex == edgeIndex) continue; // don't consider ourselves
					if (edgesVisited[nextEdgeCandidateIndex]) continue; // already part of another polygon

					// find the candidage edge dir
					ref var candidateEdge = ref straightSkeleton.edges[nextEdgeCandidateIndex];
					float2 candidateDir = math.normalize(straightSkeleton.vertices[candidateEdge.nextVertexIndex] - straightSkeleton.vertices[candidateEdge.prevVertexIndex]);

					// SignedAngle returns counter-clockwise, while we are looking for clockwise angle, hence the minus!
					float angle = -Geometry.SignedAngle(edgeDir, candidateDir);

					// we don't want to allow 180° angle switches, but any other edge closest to that is our favourite!
					if (angle > bestAngle && angle < math.PI - Geometry.EPS)
					{
						bestAngle = angle;
						nextEdgeIndex = nextEdgeCandidateIndex;
					}
				}

				if (nextEdgeIndex == -1) throw new ArgumentException($"No unused next edge found for edge starting at edge {edge}");

				return nextEdgeIndex;
			}


			[MethodImpl(MethodImplOptions.AggressiveInlining)]
			Slab CreatePolygon(int startEdgeIndex)
			{
				Slab data = new Slab();

				int edgeIndex = startEdgeIndex;
				do
				{
					ref var edge = ref straightSkeleton.edges[edgeIndex];
					data.indices.Add(edge.nextVertexIndex);
					edgeIndex = GetNextEdge(edgeIndex);
					edgesVisited[edgeIndex] = true;
				}
				while (edgeIndex != startEdgeIndex);

				return data;

			}


			// now iterate over every edge and create a polygon from it
			for (int edgeIndex = 0; edgeIndex < straightSkeleton.nEdges; ++edgeIndex)
			{
				if (edgesVisited[edgeIndex]) continue;

				polygons.Add(CreatePolygon(edgeIndex));
			}

			return polygons;
		}
	}
}