using Briganti.StraightSkeletons.Priority_Queue;
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
	public class StraightSkeletonGenerator
	{
		private WavefrontGraph wavefront;
		private VertexGraph straightSkeleton;

		private FastStructPriorityQueue<int> eventQueue;

		private float time = 0.0f;


		public StraightSkeletonGenerator(PolygonWithHoles polygonWithHoles)
		{
			// create the initial graph based on the polygon
			wavefront = new WavefrontGraph(polygonWithHoles);

			// we copy all the initial vertices & edges to the straight skeleton - they are always part
			int nStraightSkeletonVertices = wavefront.nVertices + (wavefront.nVertices - 2);
			int nStraightSkeletonEdges = (2 * wavefront.nVertices) - 3;
			straightSkeleton = new VertexGraph(nStraightSkeletonVertices, nStraightSkeletonEdges);
			for (int i = 0; i < wavefront.nVertices; ++i)
			{
				straightSkeleton.AddVertex(wavefront.vertices[i].pos, 0.0f);
			}
			for (int i = 0; i < wavefront.nEdges; ++i)
			{
				ref var edge = ref wavefront.edges[i];
				straightSkeleton.AddEdge(new Edge(edge.prevVertexIndex, edge.nextVertexIndex, edge.prevEdgeIndex, edge.nextEdgeIndex));
			}

			// we add every edge to the priority queue
			eventQueue = new FastStructPriorityQueue<int>(wavefront.maxVertices);
			for (int i = 0; i < wavefront.nEdges; ++i)
			{
				ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[i];
				eventQueue.Enqueue(i, edgeEvent.eventTime, out int queueId);
				edgeEvent.queueId = queueId;
			}

			// now we keep popping those events until we're done
			while (eventQueue.Count > 0)
			{
				int edgeIndex = eventQueue.Dequeue();
				ProcessEvent(edgeIndex);

				// TODO combine all events that happen at the same time!


			}

		}

		private void ProcessEvent(int edgeIndex)
		{
			ref Edge edge = ref wavefront.edges[edgeIndex];
			ref EdgeEvent edgeEvent = ref wavefront.edgeEvents[edgeIndex];

			// event happened somehow in the past - we skip it
			if (edgeEvent.eventTime < time) return;
			time = edgeEvent.eventTime;

			// edge event - easy
			if (edgeEvent.eventType == EventType.Edge)
			{
				ProcessEdgeEvent(edgeIndex, ref edge, ref edgeEvent);
			}
		}

		public void ProcessEdgeEvent(int eventEdgeIndex, ref Edge edge, ref EdgeEvent edgeEvent)
		{
			if (edgeEvent.eventType != EventType.Edge) throw new InvalidOperationException($"There is no edge event associated with edge {edge}.");

			// TODO first add the new vertex to the wavefront

			// TODO then add new edges to the straight skeleton


		}
	}
}