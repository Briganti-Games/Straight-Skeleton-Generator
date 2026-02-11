using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;
using Unity.Mathematics;
#if UNITY_EDITOR
using UnityEditor;
#endif
using UnityEngine;
using UnityEngine.Profiling;

namespace Briganti.StraightSkeletonGeneration
{
	public interface IEventTimeChangeListener
	{
		void AddOrUpdateEdgeEvent(int edgeIndex);
		void AddOrUpdateVertexEvent(int vertexIndex);
	}

	public class WavefrontGraph : VertexGraph
	{
		public EdgeEvent[] edgeEvents { get; private set; }
		public VertexData[] vertexDatas { get; private set; }

		private readonly IEventTimeChangeListener eventTimeListener;

		private HashSet<int> affectedEdges = new HashSet<int>();
		private HashSet<int> affectedVertices = new HashSet<int>();

		private bool debug;
		private IStraightSkeletonLogger logger;

		private float lastEventTime = 0.0f;

		private uint lastSubgraphIndex = 0;


		public WavefrontGraph(PolygonWithHoles polygonWithHoles, IEventTimeChangeListener eventTimeListener, bool debug)
		{
			this.debug = debug;
			this.eventTimeListener = eventTimeListener;

			int nVertices = polygonWithHoles.outerContour.Length;
			for (int i = 0; i < polygonWithHoles.innerContours.Count; ++i)
			{
				nVertices += polygonWithHoles.innerContours[i].Length;
			}

			int nEdges = nVertices;

			// we now add all the edges that will be added for the straight skeleton - the exact number is known
			// each vertex can spawn at most another vertex, because it either collapses an edge, or is moved when it only collapses after the max event time

			// there is a CHANCE that these values are NOT enough, because these calculations only hold for simple polygons...
			// if they contain holes, we might go over this and I couldn't find a proper formula to figure out the real max.
			int maxVertices = nVertices * 3; // we do *2 instead of *3 for no real reason - this might fix the hole issue?
			int nArcs = ((2 * nVertices) - 3);
			int maxEdges = nEdges + nArcs;

			// this can be improved and the upper limit can be calculated EXACTLY - see the lemmas in Stefan Huber's PhD
			Initialize(maxVertices, maxEdges);

			// we store additional data for each edge and vertex
			edgeEvents = new EdgeEvent[maxEdges];
			vertexDatas = new VertexData[maxVertices];
			for (int i = 0; i < maxVertices; ++i)
			{
				vertexDatas[i].queueId = -1;
				vertexDatas[i].subgraphIndex = 0;
			}

			// add all contours
			AddContour(polygonWithHoles.outerContour);
			for (int i = 0; i < polygonWithHoles.innerContours.Count; ++i)
			{
				AddContour(polygonWithHoles.innerContours[i]);
			}

			// make sure all the other edge events have their queue id reset as well
			for (int i = nEdges; i < maxEdges; ++i)
			{
				edgeEvents[i].eventType = EventType.NotInWavefront;
				edgeEvents[i].queueId = -1;
			}
		}

		public void SetLogger(IStraightSkeletonLogger logger)
		{
			this.logger = logger;
		}

		public void GenerateInitialEvents()
		{
			UpdateWavefront(0);
		}

		private void AddContour(float2[] contour)
		{
			int vertexOffset = nVertices;
			int edgeOffset = nEdges;
			for (int i = 0; i < contour.Length; ++i)
			{
				float2 curr = contour[i];

				// add the vertex
				int newVertexIndex = AddVertex(curr);

				int prevVertexIndex = vertexOffset + (i - 1 + contour.Length) % contour.Length;
				int nextVertexIndex = vertexOffset + (i + 1) % contour.Length;

				// set the vertex data
				int prevEdgeIndex = edgeOffset + (i - 1 + contour.Length) % contour.Length;
				int nextEdgeIndex = edgeOffset + i;
				vertexDatas[newVertexIndex].splitEdge = -1; // we don't split an edge by default
				vertexDatas[newVertexIndex].UpdateConnections(prevVertexIndex, nextVertexIndex, prevEdgeIndex, nextEdgeIndex);

				// also store the edge that connects this vertex and the next one
				int edgeIndex = AddEdge(new Edge(edgeOffset + i, edgeOffset + (i + 1) % contour.Length));

				// by default there is no event for this edge
				edgeEvents[edgeIndex].eventType = EventType.None;
				edgeEvents[edgeIndex].queueId = -1;

				// store the direction the edge is moving by finding the two points it connects
				float2 next = contour[(i + 1) % contour.Length];
				edgeEvents[edgeIndex].dir = math.normalize(Geometry2D.RotateMinus90Degrees(next - curr));

				affectedVertices.Add(newVertexIndex);
				affectedEdges.Add(edgeIndex);
			}
		}

		public int AddVertexToWavefrontAndRemoveEdge(int disappearedEdgeIndex, bool updateEdgeEvents)
		{
			// get the edge event
			ref EdgeEvent edgeEvent = ref edgeEvents[disappearedEdgeIndex];
			if (debug) logger?.Log($"Process edge event for {disappearedEdgeIndex}.", 3);

			// if this event has changed by now, we ignore it
			if (edgeEvent.eventType != EventType.Edge) return -1;

			float2 pos = edgeEvent.eventPos;
			float time = edgeEvent.eventTime;
			lastEventTime = Mathf.Max(lastEventTime, time);

			ref Edge disappearedEdge = ref edges[disappearedEdgeIndex];
			ref VertexData prevVertexData = ref vertexDatas[disappearedEdge.prevVertexIndex];
			ref VertexData nextVertexData = ref vertexDatas[disappearedEdge.nextVertexIndex];

			// this is a special case - we are at the end of the wavefront and it's been reduced to exactly one edge (a nook)
			// in this case, we skip this event entirely, it will be handled later
			if (prevVertexData.prevVertexIndex == disappearedEdge.nextVertexIndex && nextVertexData.nextVertexIndex == disappearedEdge.prevVertexIndex)
			{
				return -1;
			}


			if (prevVertexData.nextEdgeIndex != disappearedEdgeIndex) throw new ArgumentException($"Vertex {disappearedEdge.prevVertexIndex} is not properly linked to edge {disappearedEdgeIndex}");
			if (nextVertexData.prevEdgeIndex != disappearedEdgeIndex) throw new ArgumentException($"Vertex {disappearedEdge.nextVertexIndex} is not properly linked to edge {disappearedEdgeIndex}");

			// we flag that these vertices were removed from the wavefront, since they merged into a new vertex
			RemoveVertexFromWavefront(disappearedEdge.prevVertexIndex);
			RemoveVertexFromWavefront(disappearedEdge.nextVertexIndex);

			// the edge got permanently removed from the wavefront at this point
			RemoveEdgeFromWavefront(disappearedEdgeIndex);

			int prevEdgeIndex = prevVertexData.prevEdgeIndex;
			int nextEdgeIndex = nextVertexData.nextEdgeIndex;

			int prevVertexIndex = prevVertexData.prevVertexIndex;
			int nextVertexIndex = nextVertexData.nextVertexIndex;

			ref Edge prevEdge = ref edges[prevVertexData.prevEdgeIndex];
			ref Edge nextEdge = ref edges[nextVertexData.nextEdgeIndex];

			// firstly, just add the vertex
			int newVertexIndex = AddVertex(pos);

			// update the new vertex with the right connections
			UpdateConnections(newVertexIndex, prevVertexIndex, nextVertexIndex, prevEdgeIndex, nextEdgeIndex);
			vertexDatas[newVertexIndex].creationTime = time;

			return newVertexIndex;
		}

		public void RemoveNook(int edgeIndex)
		{
			if (debug) logger?.Log($"Remove nook edge {edgeIndex}.", 3);
			ref EdgeEvent edgeEvent = ref edgeEvents[edgeIndex];
			if (edgeEvent.eventType != EventType.Nook) throw new ArgumentException($"Edge event {edgeEvent} is not a nook event!");

			ref Edge edge = ref edges[edgeIndex];

			// the edge got permanently removed from the wavefront at this point
			RemoveEdgeFromWavefront(edgeIndex);
			lastEventTime = Mathf.Max(lastEventTime, edgeEvent.eventTime);

			// we also remove the OTHER edge of this nook from the wavefront
			ref var prevVertexData = ref vertexDatas[edge.prevVertexIndex];
			ref var otherEdgeEvent = ref edgeEvents[prevVertexData.prevEdgeIndex];
			ref var otherEdge = ref edges[prevVertexData.prevEdgeIndex];
			if (debug && (edge.prevVertexIndex != otherEdge.nextVertexIndex || edge.nextVertexIndex != otherEdge.prevVertexIndex)) throw new ArgumentException($"Edges {edge} and {otherEdge} do not form a nook!");
			RemoveEdgeFromWavefront(prevVertexData.prevEdgeIndex);

			// remove the edge from the wavefront (and along it also the other edge that forms the same nook)
			RemoveVertexFromWavefront(edge.prevVertexIndex);
			RemoveVertexFromWavefront(edge.nextVertexIndex);
		}

		public int SplitEdge(int reflexVertexIndex, int splitEdgeIndex, float time, float2 pos)
		{
			// we let the reflex vertex die - it served its purpose
			ref VertexData reflexVertexData = ref vertexDatas[reflexVertexIndex];

			// if this event has changed by now, we ignore it
			if (reflexVertexData.splitPoint != SplitPoint.Edge) throw new ArgumentException($"Split event {reflexVertexData} is not an edge split event!");

			lastEventTime = Mathf.Max(lastEventTime, time);
			if (debug) logger?.Log($"Split edge {splitEdgeIndex} at position {pos} in preparation for a split event.", 3);

			// the edge got permanently removed from the wavefront at this point
			RemoveEdgeFromWavefront(splitEdgeIndex);
			ref Edge oldEdge = ref edges[splitEdgeIndex];

			// create a new vertex at the split point
			int newVertexIndex = AddVertex(pos);
			vertexDatas[newVertexIndex].creationTime = time;

			// we need to assign a temporary velocity so we can sort this point later
			float2 prevVertex = GetVertexPosAtTime(oldEdge.prevVertexIndex, time);
			float2 nextVertex = GetVertexPosAtTime(oldEdge.nextVertexIndex, time);
			vertexDatas[newVertexIndex].velocity = math.normalize(Geometry2D.Rotate90DegreesClockwise(nextVertex - prevVertex));

			// we create two new edges from the old edge
			int newPrevEdgeIndex = AddEdge(oldEdge.prevVertexIndex, newVertexIndex);
			int newNextEdgeIndex = AddEdge(newVertexIndex, oldEdge.nextVertexIndex);

			// make sure the new edges are part of the wavefront now
			edgeEvents[newPrevEdgeIndex].eventType = EventType.None;
			edgeEvents[newNextEdgeIndex].eventType = EventType.None;
			edgeEvents[newPrevEdgeIndex].dir = edgeEvents[splitEdgeIndex].dir;
			edgeEvents[newNextEdgeIndex].dir = edgeEvents[splitEdgeIndex].dir;

			// connect the edges
			UpdateConnections(newVertexIndex, oldEdge.prevVertexIndex, oldEdge.nextVertexIndex, newPrevEdgeIndex, newNextEdgeIndex);

			// return the new vertex id
			return newVertexIndex;
		}

		public void SplitGraphAtVertices(List<int> vertexIndices, float time, float2 pos, List<Edge> newEdges, List<int> newVertexIndices)
		{
			if (vertexIndices.Count == 0) throw new ArgumentException($"Cannot do a split with 0 vertices.");
			if (debug) logger?.Log($"Split edges after vertices {string.Join(", ", vertexIndices)} collide", 3);

			int CompareReflexVertexVelocities(int vertexIndex1, int vertexIndex2)
			{
				// we sort clockwise, hence the minus before the formula
				return -Geometry2D.GetAngle(vertexDatas[vertexIndex1].velocity).CompareTo(Geometry2D.GetAngle(vertexDatas[vertexIndex2].velocity));
			}

			// sort the vertices by their velocity so that they're connected in the right order
			vertexIndices.Sort(CompareReflexVertexVelocities);
			int nVertexIndices = vertexIndices.Count;

			// First, do a prepwork step - remove any consecutive split vertices that have NO vertices in between them.
			// Since there is no new vertex loop to form between them (because there are no vertices),
			// we only keep one of these split vertices.
			for (int i = 0; i < nVertexIndices; ++i)
			{
				// all split vertices were consecutive so we don't actually continue here, since we'll be comparing the last one to itself
				if (nVertexIndices == 1) break;

				int nextListIndex = (i + 1) % nVertexIndices;
				int currVertexIndex = vertexIndices[i];
				int nextVertexIndex = vertexIndices[nextListIndex];

				ref VertexData currVertexData = ref vertexDatas[currVertexIndex];
				if (currVertexData.nextVertexIndex == nextVertexIndex)
				{

					// delete this intermediate vertex, and reconnect the chain
					RemoveVertexFromWavefront(nextVertexIndex);
					ref VertexData nextVertexData = ref vertexDatas[nextVertexIndex];
					UpdateConnections(currVertexIndex, currVertexData.prevVertexIndex, nextVertexData.nextVertexIndex, currVertexData.prevEdgeIndex, nextVertexData.nextEdgeIndex);

					// remove from the list of split vertices
					// the fact that we might remove index 0 is not a real issue 
					vertexIndices.RemoveAt((i + 1) % nVertexIndices);
					--i;
					if (nextListIndex == 0) --i; // we removed one ahead of us so we do another step back
					--nVertexIndices;
				}
			}

			// we're going to remove all these vertices and replace them by new ones that connect the parts between them
			for (int i = 0; i < nVertexIndices; ++i)
			{
				// get the original vertex data
				int currVertexIndex = vertexIndices[i];
				ref VertexData currVertexData = ref vertexDatas[currVertexIndex];

				// we kill this vertex and spawn a new one at the event position
				RemoveVertexFromWavefront(currVertexIndex);

				// spawn a new vertex at the event position
				int newVertexIndex = AddVertex(pos);
				newVertexIndices.Add(newVertexIndex);
				vertexDatas[newVertexIndex].creationTime = time;

				// if the old vertex is a reflex vertex or we are a endpoint split, we add a new edge
				// if were an edge split, we just created a vertex on the existing edge and we don't want to
				// draw an edge from that vertex to the new one at the exact same pos.
				float2 prevPos = vertices[currVertexIndex];
				float2 newPos = vertices[newVertexIndex];
				if (math.distancesq(prevPos, newPos) > Geometry2D.EPSSQ)
				{
					newEdges.Add(new Edge(currVertexIndex, newVertexIndex));
				}

				// we now link the new vertex to the previous and next connections of the previous and next vertices involved in this business
				//ref VertexData prevVertexData = ref vertexDatas[vertexIndices[(i-1+nVertexIndices)%nVertexIndices]];
				int nextVertexIndex = vertexIndices[(i + 1) % nVertexIndices];
				ref VertexData nextVertexData = ref vertexDatas[nextVertexIndex];
				UpdateConnections(newVertexIndex, nextVertexData.prevVertexIndex, currVertexData.nextVertexIndex, nextVertexData.prevEdgeIndex, currVertexData.nextEdgeIndex);
			}

			// now go over all new vertex indices, and assign a new subgraph to the loop that this index is part of - separate subgraphs are not allowed to interact with each other!
			// We skip the first one, because we can keep that one's subgraph, since we're changing all the other ones. Just a minor optimization.
			for (int i = 1; i < newVertexIndices.Count; ++i)
			{
				AssignNewSubgraphToLoop(newVertexIndices[i]);
			}
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private void UpdateConnections(int vertexIndex, int prevVertexIndex, int nextVertexIndex, int prevEdgeIndex, int nextEdgeIndex)
		{
			// update the vertex itself
			ref var vertexData = ref vertexDatas[vertexIndex];

			// the old connected edges are affected since they might be disconnected
			if (vertexData.type != WavefrontVertexType.Unknown)
			{
				affectedEdges.Add(vertexData.prevEdgeIndex);
				affectedEdges.Add(vertexData.nextEdgeIndex);
			}

			vertexData.UpdateConnections(prevVertexIndex, nextVertexIndex, prevEdgeIndex, nextEdgeIndex);

			// also update the adjacent edges
			edges[vertexData.prevEdgeIndex].nextVertexIndex = vertexIndex;
			edges[vertexData.nextEdgeIndex].prevVertexIndex = vertexIndex;

			// copy the subgraph index of an adjacent vertex
			vertexData.subgraphIndex = vertexDatas[prevVertexIndex].subgraphIndex;

			// also update the adjacent vertices
			vertexDatas[prevVertexIndex].nextVertexIndex = vertexIndex;
			vertexDatas[nextVertexIndex].prevVertexIndex = vertexIndex;
			vertexDatas[prevVertexIndex].nextEdgeIndex = prevEdgeIndex;
			vertexDatas[nextVertexIndex].prevEdgeIndex = nextEdgeIndex;

			// remember these are affected by graph changes, so we can update them later
			affectedVertices.Add(vertexIndex);

			// the newly connected edges are affected
			affectedEdges.Add(prevEdgeIndex);
			affectedEdges.Add(nextEdgeIndex);
		}

		private void AssignNewSubgraphToLoop(int startVertexIndex)
		{
			++lastSubgraphIndex;

			// if we somehow reached the very last subgraph, we reset everything and assign new indices to the entire graph
			if (lastSubgraphIndex == uint.MaxValue)
			{
				lastSubgraphIndex = 0;

				for (int i = 0; i < nVertices; ++i)
				{
					vertexDatas[i].subgraphIndex = 0;
				}

				for (int i = 0; i < nVertices; ++i)
				{
					if (vertexDatas[i].subgraphIndex == 0) AssignNewSubgraphToLoop(i);
				}

				return;
			}

			// otherwise, just assign to this loop
			int vertexIndex = startVertexIndex;
			do
			{
				vertexDatas[vertexIndex].subgraphIndex = lastSubgraphIndex;
				vertexIndex = vertexDatas[vertexIndex].nextVertexIndex;
			}
			while (vertexIndex != startVertexIndex);
			
		}

		private void RemoveVertexFromWavefront(int vertexIndex)
		{
			ref VertexData vertexData = ref vertexDatas[vertexIndex];
			vertexData.inWavefront = false;
			vertexData.partOfSplitEvent = false;
			eventTimeListener.AddOrUpdateVertexEvent(vertexIndex);
		}

		private void RemoveEdgeFromWavefront(int edgeIndex)
		{
			ref EdgeEvent edgeEvent = ref edgeEvents[edgeIndex];
			edgeEvent.eventType = EventType.NotInWavefront;
		}


		public void UpdateWavefront(float time)
		{
			/*UpdateEntireWavefront(time);
			return;*/

			// first, make sure all vertices know their velocity etc
			foreach (int vertexIndex in affectedVertices)
			{
				UpdateVertexData(vertexIndex, time);
			}

			// we update the edge/nook events for all edges
			foreach (int edgeIndex in affectedEdges)
			{
				UpdateEdgeEventTime(edgeIndex, time);
			}

			// now for each reflex vertex, calculate the edge they will be splitting and update that edge's event if that's the first thing that'll happen to it
			bool changed;
			for (int vertexIndex = 0; vertexIndex < nVertices; ++vertexIndex)
			{
				changed = false;
				ref var vertexData = ref vertexDatas[vertexIndex];


				if (vertexData.inWavefront && vertexData.type == WavefrontVertexType.Reflex)
				{
					if (vertexData.partOfSplitEvent)
					{

						// if our current split event is for a deleted edge, we forget it!
						if (affectedEdges.Contains(vertexData.splitEdge) || edgeEvents[vertexData.splitEdge].eventType == EventType.NotInWavefront)
						{
							vertexData.partOfSplitEvent = false;
							changed = true;
						}
					}

					// if this is an affected vertex, we need to recheck all edges
					if (affectedVertices.Contains(vertexIndex))
					{
						// always changed because we reset our state
						changed = true;

						for (int edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex)
						{
							if (edgeEvents[edgeIndex].eventType != EventType.NotInWavefront)
							{
								Profiler.BeginSample("UpdateSplitEventTime");
								if (UpdateSplitEventTime(vertexIndex, edgeIndex, time)) changed = true;
								Profiler.EndSample();
							}
						}
					}

					// if the vertex itself did not change, we do need to check all affected edges
					else
					{
						foreach (int edgeIndex in affectedEdges)
						{
							if (edgeEvents[edgeIndex].eventType != EventType.NotInWavefront)
							{
								Profiler.BeginSample("UpdateSplitEventTime");
								if (UpdateSplitEventTime(vertexIndex, edgeIndex, time)) changed = true;
								Profiler.EndSample();
							}
						}
					}

					// now assign the split event to the edge event if there is one
					if (changed)
					{
						eventTimeListener.AddOrUpdateVertexEvent(vertexIndex);
					}
				}
			}

			affectedVertices.Clear();
			affectedEdges.Clear();
		}

		public void UpdateEntireWavefront(float time)
		{
			for (int i = 0; i < nVertices; ++i)
			{
				UpdateVertexData(i, time);
			}

			for (int i = 0; i < nEdges; ++i)
			{
				UpdateEdgeEventTime(i, time);
			}

			// now for each reflex vertex, calculate the edge they will be splitting and update that edge's event if that's the first thing that'll happen to it
			for (int vertexIndex = 0; vertexIndex < nVertices; ++vertexIndex)
			{
				if (vertexDatas[vertexIndex].inWavefront && vertexDatas[vertexIndex].type == WavefrontVertexType.Reflex)
				{
					vertexDatas[vertexIndex].partOfSplitEvent = false;

					for (int edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex)
					{
						if (edgeEvents[edgeIndex].eventType != EventType.NotInWavefront)
						{
							Profiler.BeginSample("UpdateSplitEventTime");
							UpdateSplitEventTime(vertexIndex, edgeIndex, time);
							Profiler.EndSample();
						}
					}

					// now assign the split event to the edge event if there is one
					eventTimeListener.AddOrUpdateVertexEvent(vertexIndex);
				}
			}

			affectedVertices.Clear();
			affectedEdges.Clear();
		}

		public int GetNextVertexIndex(int vertexIndex)
		{
			ref VertexData vertexData = ref vertexDatas[vertexIndex];
			ref Edge edge = ref edges[vertexData.nextEdgeIndex];
			return edge.nextVertexIndex;
		}

		public int GetPrevVertexIndex(int vertexIndex)
		{
			ref VertexData vertexData = ref vertexDatas[vertexIndex];
			ref Edge edge = ref edges[vertexData.prevEdgeIndex];
			return edge.prevVertexIndex;
		}

		private void UpdateVertexData(int vertexIndex, float time)
		{
			Profiler.BeginSample("UpdateVertexData");

			ref VertexData vertexData = ref vertexDatas[vertexIndex];

			// we reset our split event state - and also the edge that was associated with our split
			// because of a different event, this vertex might have been removed but might have still been involved in a split, so we need to remove that one
			if (vertexData.partOfSplitEvent)
			{
				vertexData.partOfSplitEvent = false;
			}

			// if we're not part of the wavefront we don't care anymore about the rest
			if (!vertexData.inWavefront)
			{
				Profiler.EndSample();
				return;
			}

			// we never re-calculate a vertex velocity (is this right?)
			if (vertexData.type == WavefrontVertexType.Unknown)
			{
				float2 prevToCurrLineDir = edgeEvents[vertexData.prevEdgeIndex].dir;
				float2 currToNextLineDir = edgeEvents[vertexData.nextEdgeIndex].dir;

				float2 velocity = CalculateVelocity(prevToCurrLineDir, currToNextLineDir);
				WavefrontVertexType type = WavefrontVertexType.Convex;
				if (math.length(velocity) > Geometry2D.EPS)
				{
					float angleFromPrevToNextLine = Geometry2D.SignedAngle(prevToCurrLineDir, currToNextLineDir);
					bool isReflex = angleFromPrevToNextLine >= -Geometry2D.EPS;

					type = (isReflex ? WavefrontVertexType.Reflex : WavefrontVertexType.Convex);
				}

				vertexData.velocity = velocity;
				vertexData.type = type;
			}

			Profiler.EndSample();
		}

		public static float2 CalculateVelocity(in float2 prevToCurrLineDir, in float2 currToNextLineDir)
		{
			// if the adjacent move dirs cancel each other out, we don't move
			float2 sum = prevToCurrLineDir + currToNextLineDir;
			if (math.lengthsq(sum) < Geometry2D.EPSSQ)
			{
				return float2.zero;
			}

			float2 dir = math.normalize(prevToCurrLineDir + currToNextLineDir);
			float speed = 2.0f / (math.dot(prevToCurrLineDir, dir) + math.dot(currToNextLineDir, dir));

			float2 velocity = dir * speed;
			//Debug.Log("Vertex #" + vertexIndex + " has dir " + dir + " and speed " + speed + " based on " + moveDir1 + " and " + moveDir2 + " and final dir has angle " + (Geometry.GetAngle(dir) * Mathf.Rad2Deg) + ", angle between two adjacent dirs is " + (Geometry.SignedAngle(prevToCurr, nextToCurr) * Mathf.Rad2Deg));
			return velocity;
		}

		public void UpdateEdgeEventTime(int edgeIndex, float time)
		{
			// default to max
			ref Edge edge = ref edges[edgeIndex];
			ref EdgeEvent eventData = ref edgeEvents[edgeIndex];

			// once we're removed from the wavefront, nothing happens to this edge anymore
			if (eventData.eventType == EventType.NotInWavefront)
			{
				// remove from the eventqueue
				eventTimeListener.AddOrUpdateEdgeEvent(edgeIndex);
				return;
			}

			eventData.Reset();

			// the edge must be either 100% in the current wavefront or removed from it
			bool vertexInWavefront1 = vertexDatas[edge.prevVertexIndex].inWavefront;
			bool vertexInWavefront2 = vertexDatas[edge.nextVertexIndex].inWavefront;
			if (vertexInWavefront1 != vertexInWavefront2) throw new InvalidOperationException($"Edge {edgeIndex} {edge} is partly in the wavefront and partly not, because it connects vertices {vertices[edge.prevVertexIndex]} {vertexDatas[edge.prevVertexIndex]} and {vertices[edge.nextVertexIndex]} {vertexDatas[edge.nextVertexIndex]}. This is illegal!");
			if (!vertexInWavefront1)
			{
				eventTimeListener.AddOrUpdateEdgeEvent(edgeIndex);
				return;
			}

			// first, see if we're at the end if a wavefront and it has collapsed to one edge - a nook event
			bool wasNook = UpdateNookEventTime(ref edge, ref eventData, time);
			if (!wasNook)
			{

				// then, the easy part - calculate the edge event
				UpdateEdgeEventTime(ref edge, ref eventData, time);
			}

			// add to the queue
			eventTimeListener.AddOrUpdateEdgeEvent(edgeIndex);
		}

		private bool UpdateNookEventTime(ref Edge edge, ref EdgeEvent eventData, float time)
		{
			Profiler.BeginSample("UpdateNookEventTime");
			ref VertexData prevVertexData = ref vertexDatas[edge.prevVertexIndex];
			ref VertexData nextVertexData = ref vertexDatas[edge.nextVertexIndex];

			// this is a special case - we are at the end of the wavefront and it's been reduced to exactly one edge (a nook)
			// in this case, we skip this event entirely, it will be handled later
			if (prevVertexData.prevVertexIndex == edge.nextVertexIndex && nextVertexData.nextVertexIndex == edge.prevVertexIndex)
			{
				// because a nook event always comes in pairs, we check whether the OTHER edge has already been flagged as a nook
				// in that wase, we just skip this one
				ref var otherEdgeEvent = ref edgeEvents[prevVertexData.prevEdgeIndex];
				if (otherEdgeEvent.eventType == EventType.Nook)
				{
					eventData.eventType = EventType.None;
					eventData.eventTime = float.MaxValue;
					eventData.eventPos = float2.zero;
				}

				// first of the two edges we encountered - it's a nook!
				else
				{
					eventData.eventType = EventType.Nook;
					eventData.eventTime = time;
					eventData.eventPos = float2.zero;
				}

				Profiler.EndSample();
				return true;
			}

			Profiler.EndSample();
			return false;
		}

		private void UpdateEdgeEventTime(ref Edge edge, ref EdgeEvent eventData, float _)
		{
			Profiler.BeginSample("UpdateEdgeEventTime");

			// we store the edge event time in the prevVertex associated with the edge!
			float2 prevVertex = vertices[edge.prevVertexIndex];
			float2 nextVertex = vertices[edge.nextVertexIndex];
			ref VertexData prevData = ref vertexDatas[edge.prevVertexIndex];
			ref VertexData nextData = ref vertexDatas[edge.nextVertexIndex];

			[MethodImpl(MethodImplOptions.AggressiveInlining)]
			bool IsCloseEnoughToPositive(float t0, float2 velocity0, float t1, float2 velocity1)
			{
				// if this SHOULD be false, we calculate how far we actually travelled, and if it's not that far (because velocity is very small), we let it slip
				if (t0 <= -Geometry2D.EPS)
				{
					float2 offset0 = velocity0 * t0;
					float d0 = math.lengthsq(offset0);
					if (d0 > Geometry2D.EPSSQ_LOWPRECISION) return false;
				}
				if (t1 <= -Geometry2D.EPS)
				{
					float2 offset1 = velocity1 * t1;
					float d1 = math.lengthsq(offset1);
					if (d1 > Geometry2D.EPSSQ_LOWPRECISION) return false;
				}

				return true;
			}

			// this is a very simple case - both vertices are at the exact same position, so the edge is basically already zero
			// this can be caused by two edge events happening at the same position, and causing two new vertices there.
			if (math.distancesq(prevVertex, nextVertex) < Geometry2D.EPS)
			{
				eventData.eventType = EventType.Edge;
				eventData.eventPos = prevVertex;
				eventData.eventTime = prevData.creationTime;
				Profiler.EndSample();
				return;
			}

			if (Geometry2D.GetLineIntersection(prevVertex, prevVertex + prevData.velocity, nextVertex, nextVertex + nextData.velocity, out float t0, out float t1))
			{
				if (IsCloseEnoughToPositive(t0, prevData.velocity, t1, nextData.velocity))
				{
					eventData.eventType = EventType.Edge;
					eventData.eventPos = prevVertex + prevData.velocity * t0;

					float time0 = prevData.creationTime + t0;
					float time1 = nextData.creationTime + t1;

					float eventTime = Mathf.Max(time0, time1);
					eventData.eventTime = eventTime;
				}
			}

			Profiler.EndSample();
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		public float2 GetVertexPosAtTime(int vertexIndex, float time)
		{
			float2 vertex = vertices[vertexIndex];
			ref VertexData vertexData = ref vertexDatas[vertexIndex];
			return GetVertexPosAtTime(vertex, vertexData, time);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private float2 GetVertexPosAtTime(in float2 vertex, in VertexData vertexData, float time)
		{
			//if (time < vertexData.creationTime - Geometry.EPS) throw new ArgumentException($"Time {time} lies in the past of the creation point of vertex {vertex} ({vertexData.creationTime}), this should not be possible.");
			if (time < vertexData.creationTime - Geometry2D.EPS) time = vertexData.creationTime;
			time -= vertexData.creationTime;
			return vertex + vertexData.velocity * time;
		}

		private bool UpdateSplitEventTime(int reflexVertexIndex, int edgeIndex, float currentTime)
		{
			ref Edge edge = ref edges[edgeIndex];

			// can't be adjacent
			if (edge.prevVertexIndex == reflexVertexIndex || edge.nextVertexIndex == reflexVertexIndex) return false;
			if (edgeEvents[edgeIndex].eventType == EventType.NotInWavefront) return false;

			float2 prevVertex = GetVertexPosAtTime(edge.prevVertexIndex, currentTime);
			float2 nextVertex = GetVertexPosAtTime(edge.nextVertexIndex, currentTime);

			float2 reflexVertex = GetVertexPosAtTime(reflexVertexIndex, currentTime);
			ref VertexData vertexData = ref vertexDatas[reflexVertexIndex];
			if (!vertexData.inWavefront) return false;
			if (vertexData.subgraphIndex != GetEdgeSubgraphIndex(edgeIndex)) return false; // we are in a different subgraph of the wavefront, caused by a split!
			if (vertexData.type != WavefrontVertexType.Reflex) throw new ArgumentException($"Vertex {reflexVertexIndex} is not a reflex vertex.");
			float2 reflexVelocity = vertexData.velocity;

			// we calculate the collision point between the reflex vertex and the edge
			bool collisionFound = CalculateCollisionPoint(reflexVertex, reflexVelocity, prevVertex, nextVertex, out float L);

			// no collision found
			if (!collisionFound) return false;

			// we can now determine the position where the edge meets the reflex vertex
			//float distanceAlongVertexLine = L * V;
			float time = L;
			float2 eventPos = reflexVertex + reflexVelocity * time;

			// because lines always move at constant speed, the time is identical to the move distance of the edge
			float eventTime = currentTime + time;

			// if, at the exact same time, the reflex vertex is reduced because it is also part of an edge event, we ignore this split,
			// because it will be invalid, since the reflex vertex will be desolved in the split event before it can split anything.
			if (IsDestroyedByEdgeEventBefore(vertexData, eventTime)) return false;

			// if the reflex vertex lies ON the prev or next vertex, we have an immediate-time vertex split
			if (math.distancesq(prevVertex, reflexVertex) < Geometry2D.EPSSQ || math.distancesq(nextVertex, reflexVertex) < Geometry2D.EPSSQ)
			{
				// if the vertex SPAWNED on the line and immediately triggers a collision,
				// it means we previously already processed this collision. We ignore it!
				float2 originalVertexPosition = vertices[reflexVertexIndex];
				if (math.distancesq(prevVertex, originalVertexPosition) < Geometry2D.EPSSQ || math.distancesq(nextVertex, originalVertexPosition) < Geometry2D.EPSSQ)
				{
					return false;
				}

				// this is the best split (earliest and preferably edge before prev/next) for this vertex, so we assign it!
				var prevNextSplitPoint = (math.distancesq(prevVertex, reflexVertex) < Geometry2D.EPSSQ ? SplitPoint.PrevVertex : SplitPoint.NextVertex);
				AssignSplitToReflexVertex(edgeIndex, ref vertexData, currentTime, eventPos, prevNextSplitPoint);
				return true;
			}

			// we had a split for this reflex vertex earlier - this one is irrelevant
			if (vertexData.partOfSplitEvent && eventTime >= vertexData.splitTime + Geometry2D.EPS) return false;

			// this edge already has an event earlier
			// if there is one at the same time, we want to overwrite it with a split event if it exists
			ref EdgeEvent edgeEvent = ref edgeEvents[edgeIndex];
			if (edgeEvent.eventTime < eventTime - Geometry2D.EPS) return false;

			// now finally, we need to make sure that this still falls into the edge once it moved that far
			float2 prevVertexAtEventTime = GetVertexPosAtTime(edge.prevVertexIndex, eventTime);
			float2 nextVertexAtEventTime = GetVertexPosAtTime(edge.nextVertexIndex, eventTime);
			Geometry2D.ProjectPointOnLine(prevVertexAtEventTime, nextVertexAtEventTime, eventPos, out float tLine);
			if (tLine < -Geometry2D.EPS || tLine > 1 + Geometry2D.EPS || tLine == float.NaN) return false;

			// if we split at exactly the corner, we have a special type of split event - this event will NOT generate new edges,
			// but will still split the entire graph into two pieces and spawn new vertices at the split point.
			int splitPointVertexIndex = -1;
			SplitPoint splitPoint = SplitPoint.Edge;
			if (tLine < Geometry2D.EPS)
			{
				splitPointVertexIndex = edge.prevVertexIndex;
				splitPoint = SplitPoint.PrevVertex;
			}
			else if (tLine > 1.0f - Geometry2D.EPS)
			{
				splitPointVertexIndex = edge.nextVertexIndex;
				splitPoint = SplitPoint.NextVertex;
			}

			// we prioritize edge splits and convex vertex splits over others
			if (!IsBetterSplit(vertexData, eventTime, splitPoint, splitPointVertexIndex)) return false;

			AssignSplitToReflexVertex(edgeIndex, ref vertexData, eventTime, eventPos, splitPoint);
			return true;
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private uint GetEdgeSubgraphIndex(int edgeIndex)
		{
			int prevVertexIndex = edges[edgeIndex].prevVertexIndex;
			ref VertexData vertexData = ref vertexDatas[prevVertexIndex];
			return vertexData.subgraphIndex;
		}

		private void AssignSplitToReflexVertex(int edgeIndex, ref VertexData vertexData, float eventTime, float2 eventPos, SplitPoint splitPoint)
		{
			vertexData.partOfSplitEvent = true;
			vertexData.splitPos = eventPos;
			vertexData.splitTime = eventTime;
			vertexData.splitEdge = edgeIndex;
			vertexData.splitPoint = splitPoint;
		}

		private bool IsDestroyedByEdgeEventBefore(in VertexData vertexData, float eventTime)
		{
			ref EdgeEvent prevEdge = ref edgeEvents[vertexData.prevEdgeIndex];
			if (prevEdge.eventType == EventType.Edge && prevEdge.eventTime < eventTime + Geometry2D.EPS) return true;

			ref EdgeEvent nextEdge = ref edgeEvents[vertexData.nextEdgeIndex];
			if (nextEdge.eventType == EventType.Edge && nextEdge.eventTime < eventTime + Geometry2D.EPS) return true;

			return false;
		}

		private bool IsBetterSplit(in VertexData currentSplit, float newSplitEventTime, SplitPoint newSplitPoint, int newSplitPointVertexIndex)
		{

			// if there is no current split, it is always better
			if (!currentSplit.partOfSplitEvent) return true;

			// if we are earlier than the previous split, we definitely are better
			if (newSplitEventTime < currentSplit.splitTime - Geometry2D.EPS) return true;

			// if the current split is NOT an edge split and the new one is, the new one is better
			if (currentSplit.splitPoint != SplitPoint.Edge && newSplitPoint == SplitPoint.Edge) return true;

			// if none of the split points are edges, we prefer the convex vertices (which represent non-reflex edges we're going to collide with)
			if (currentSplit.splitPoint != SplitPoint.Edge && newSplitPoint != SplitPoint.Edge)
			{

				ref VertexData oldSplitPointData = ref vertexDatas[edges[currentSplit.splitEdge].GetPoint(currentSplit.splitPoint)];
				ref VertexData newSplitPointData = ref vertexDatas[newSplitPointVertexIndex];
				if (oldSplitPointData.type != WavefrontVertexType.Convex && newSplitPointData.type == WavefrontVertexType.Convex) return true;
			}

			return false;
		}

		private bool CalculateCollisionPoint(float2 p, float2 v, float2 p1, float2 p2, out float y)
		{
			y = 0.0f;

			// if the two points are identical, we don't have a line we can properly calculate stuff on!
			if (math.distancesq(p1, p2) < Geometry2D.EPSSQ) return false;

			float2 normal = Geometry2D.Rotate90DegreesClockwise(p2 - p1);

			// first, we offset everthing so that p1 is the origin, ensuring that the line goes through the origin
			p -= p1;

			// now we need to rotate everything so that the normal is pointing up
			float angle = Geometry2D.SignedAngle(normal, new float2(0, 1));

			// rotate the reflex velocity
			v = Geometry2D.Rotate(v, angle);
			p = Geometry2D.Rotate(p, angle);

			// if we are at the line, we meet immediately, no matter our speed
			if (math.abs(p.y) < Geometry2D.EPS) return true;

			// if we are moving at exactly the speed of the line, but we are not at the line at the start, we also never meet because we'll be moving parallel
			if (math.abs(v.y - 1) < Geometry2D.EPS) return false;

			// if we are below the line and not moving upwards fast enough, we never meet
			if (p.y < -Geometry2D.EPS && v.y <= 1 + Geometry2D.EPS) return false;

			// if we are above the line and not moving SLOW enough, we never meet
			if (p.y > Geometry2D.EPS && v.y >= 1 - Geometry2D.EPS) return false;

			// now calculate the actual moving point
			y = p.y / (1 - v.y);
			return true;
		}

		public float GetVertexTime(int vertexIndex)
		{
			return vertexDatas[vertexIndex].creationTime;
		}

		public bool IsPartOfWavefront(int edge)
		{
			return edgeEvents[edge].eventType != EventType.NotInWavefront;
		}

		public void ValidateState(float time)
		{
			if (!debug) return;

			for (int i = 0; i < nVertices; ++i)
			{
				ref var vertexData = ref vertexDatas[i];
				if (!vertexData.inWavefront) continue;

				ref var prevEdge = ref edges[vertexData.prevEdgeIndex];
				float2 prevEdgeDir = edgeEvents[vertexData.prevEdgeIndex].dir;
				ref var nextEdge = ref edges[vertexData.nextEdgeIndex];
				float2 nextEdgeDir = edgeEvents[vertexData.nextEdgeIndex].dir;

				if (prevEdge.nextVertexIndex != i) throw new ArgumentException($"Prev edge {prevEdge} (idx {vertexData.prevEdgeIndex}) does not connect correctly with vertex {vertexData} (idx {i})");
				if (nextEdge.prevVertexIndex != i) throw new ArgumentException($"Next edge {nextEdge} (idx {vertexData.nextEdgeIndex}) does not connect correctly with vertex {vertexData} (idx {i})");
				if (prevEdge.prevVertexIndex != vertexData.prevVertexIndex) throw new ArgumentException($"Prev edge {prevEdge} (idx {vertexData.prevEdgeIndex}) does not have the same previous vertex {prevEdge.prevVertexIndex} defined as the vertex {vertexData} (idx {i})");
				if (nextEdge.nextVertexIndex != vertexData.nextVertexIndex) throw new ArgumentException($"Nexts edge {nextEdge} (idx {vertexData.nextEdgeIndex}) does not have the same previous vertex {prevEdge.nextVertexIndex} defined as the vertex {vertexData} (idx {i})");

				// make sure the vertex moves in the same direction as both adjacent edges - we can never be moving in opposite directions!
				float2 velocity = vertexData.velocity;
				if (math.dot(velocity, prevEdgeDir) < -Geometry2D.EPS) throw new ArgumentException($"Prev edge {prevEdge} has move dir {prevEdgeDir}, but adjacent vertex {i} is moving in opposite direction {velocity}.");
				if (math.dot(velocity, nextEdgeDir) < -Geometry2D.EPS) throw new ArgumentException($"Next edge {nextEdge} has move dir {nextEdgeDir}, but adjacent vertex {i} is moving in opposite direction {velocity}.");
			}

			PrintState(time);

		}

		private void PrintState(float time)
		{
			string GetVertexDescription(int vertexIndex)
			{
				return $"{vertexIndex} ({GetVertexPosAtTime(vertexIndex, time)})";
			}

			// print the vertex chains
			StringBuilder ss = new StringBuilder();
			ss.AppendLine("Chains at time " + time + ":");
			int chainCount = 0;
			HashSet<int> visitedVertices = new HashSet<int>();
			for (int i = 0; i < nVertices; ++i)
			{
				ref var vertexData = ref vertexDatas[i];
				if (!vertexData.inWavefront) continue;
				if (visitedVertices.Contains(i)) continue;

				int initialVertexIndex = i;
				int vertexIndex = initialVertexIndex;
				ss.Append("Chain #" + (++chainCount) + ":");
				ss.Append(GetVertexDescription(initialVertexIndex));
				visitedVertices.Add(initialVertexIndex);
				int nVerticesVisited = 0;
				do
				{
					vertexIndex = vertexData.nextVertexIndex;
					vertexData = ref vertexDatas[vertexIndex];
					ss.Append(" -> ");
					ss.Append(GetVertexDescription(vertexIndex));
					visitedVertices.Add(vertexIndex);
					++nVerticesVisited;
				}
				while (vertexIndex != initialVertexIndex && nVerticesVisited < nVertices);


				if (nVerticesVisited >= nVertices)
				{
					ss.Append(" -> INFINITE LOOP");
					throw new ArgumentException($"There is a chain in the wavefront that is NOT a loop, this is an invalid state:\n" + ss.ToString());
				}

				ss.AppendLine();
			}


			logger?.Log(ToString(), 0);
		}

		public override string ToString()
		{
			float time = lastEventTime;

			string GetVertexDescription(int vertexIndex)
			{
				return $"{vertexIndex}";
			}

			// print the vertex chains
			StringBuilder ss = new StringBuilder();
			int chainCount = 0;
			HashSet<int> visitedVertices = new HashSet<int>();
			for (int i = 0; i < nVertices; ++i)
			{
				ref var vertexData = ref vertexDatas[i];
				if (!vertexData.inWavefront) continue;
				if (visitedVertices.Contains(i)) continue;

				int initialVertexIndex = i;
				int vertexIndex = initialVertexIndex;
				ss.Append("#" + (++chainCount) + ":");
				ss.Append(GetVertexDescription(initialVertexIndex));
				visitedVertices.Add(initialVertexIndex);
				int nVerticesVisited = 0;
				do
				{
					vertexIndex = vertexData.nextVertexIndex;
					vertexData = ref vertexDatas[vertexIndex];
					ss.Append(" ");
					ss.Append(GetVertexDescription(vertexIndex));
					visitedVertices.Add(vertexIndex);
					++nVerticesVisited;
				}
				while (vertexIndex != initialVertexIndex && nVerticesVisited < nVertices);


				if (nVerticesVisited >= nVertices)
				{
					ss.Append(" -> INFINITE LOOP");
					throw new ArgumentException($"There is a chain in the wavefront that is NOT a loop, this is an invalid state:\n" + ss.ToString());
				}

				ss.Append(" ; ");
			}

			return ss.ToString();
		}

#if UNITY_EDITOR
		public static GUIStyle vertexGuiStyle = new GUIStyle() { normal = new GUIStyleState() { textColor = Color.cyan }, fontSize = 18, alignment = TextAnchor.MiddleCenter, fontStyle = FontStyle.Bold };
		public static GUIStyle edgeGuiStyle = new GUIStyle() { normal = new GUIStyleState() { textColor = Color.red }, fontSize = 18, alignment = TextAnchor.MiddleCenter, fontStyle = FontStyle.Bold };
#endif

		public void DrawGizmos(float time, Camera camera)
		{
			Gizmos.color = Color.cyan;

			// calculate the distance from the camera of all the points, so we can apply a scale for some offsets based on the distance from the camera
			float minDistanceFromCamera = float.MaxValue;
			Vector3 cameraPos = Gizmos.matrix.MultiplyPoint(camera.transform.position);
			for (int i = 0; i < nVertices; ++i)
			{
				ref var vertexData = ref vertexDatas[i];
				var vertex = GetVertexPosAtTime(i, vertexData.creationTime);
				minDistanceFromCamera = Mathf.Min(minDistanceFromCamera, Vector3.Distance(new Vector3(vertex.x, vertexData.creationTime, vertex.y), cameraPos));
			}

			// based on the distance from the camera, we can now calculate an appropriate offset for all the labels
			float viewportOffset = 0.005f;
			Ray ray1 = camera.ViewportPointToRay(new Vector2(0, 0));
			Ray ray2 = camera.ViewportPointToRay(new Vector2(viewportOffset, 0));
			Vector3 p1 = ray1.GetPoint(minDistanceFromCamera);
			Vector3 p2 = ray2.GetPoint(minDistanceFromCamera);
			float offset = Vector3.Distance(p1, p2);

			for (int i = 0; i < nEdges; ++i)
			{
				ref var edge = ref edges[i];
				if (vertexDatas[edge.prevVertexIndex].inWavefront && vertexDatas[edge.prevVertexIndex].nextEdgeIndex == i && vertexDatas[edge.nextVertexIndex].inWavefront && vertexDatas[edge.nextVertexIndex].prevEdgeIndex == i)
				{
					float prevCreationTime = vertexDatas[edge.prevVertexIndex].creationTime;
					float nextCreationTime = vertexDatas[edge.nextVertexIndex].creationTime;

					prevCreationTime = time;
					nextCreationTime = time;

					var prevVertex = GetVertexPosAtTime(edge.prevVertexIndex, prevCreationTime);
					var nextVertex = GetVertexPosAtTime(edge.nextVertexIndex, nextCreationTime);

					Gizmos.DrawLine(new Vector3(prevVertex.x, prevCreationTime, prevVertex.y), new Vector3(nextVertex.x, nextCreationTime, nextVertex.y));

#if UNITY_EDITOR
					float2 midPoint = (prevVertex + nextVertex) * 0.5f;
					float2 normalBack = math.normalize(Geometry2D.Rotate90DegreesCounterClockwise(nextVertex - prevVertex));
					midPoint += normalBack * offset; // 0.1f
					float midCreationTime = (prevCreationTime + nextCreationTime) * 0.5f;
					Handles.Label(Gizmos.matrix.MultiplyPoint(new Vector3(midPoint.x, midCreationTime, midPoint.y)), "" + i, edgeGuiStyle);
#endif
				}
			}

			Gizmos.color = new Color(0.5f, 0.5f, 1.0f, 1f);

			for (int i = 0; i < nVertices; ++i)
			{
				ref var vertexData = ref vertexDatas[i];
				//var vertex = GetVertexPosAtTime(i, time);
				var vertex = GetVertexPosAtTime(i, vertexData.creationTime);
				if (!vertexData.inWavefront) continue;

				var endPoint = vertex + math.normalize(vertexData.velocity) * offset * 3;
				float creationTime = vertexDatas[i].creationTime;

				Gizmos.DrawLine(new Vector3(vertex.x, creationTime, vertex.y), new Vector3(endPoint.x, creationTime, endPoint.y));

#if UNITY_EDITOR
				var idPoint = vertex;
				if (math.length(vertexData.velocity) > Geometry2D.EPS) idPoint -= math.normalize(vertexData.velocity) * offset;
				Handles.Label(Gizmos.matrix.MultiplyPoint(new Vector3(idPoint.x, creationTime, idPoint.y)), "" + i, vertexGuiStyle);
#endif
			}
		}

		protected override void IncreaseVertexCapacity(int nVertices)
		{
			base.IncreaseVertexCapacity(nVertices);

			int oldLength = vertexDatas.Length;
			VertexData[] newVertexDatas = new VertexData[nVertices];
			Array.Copy(vertexDatas, newVertexDatas, oldLength);
			vertexDatas = newVertexDatas;

			for (int i = oldLength; i < maxVertices; ++i)
			{
				vertexDatas[i].queueId = -1;
			}
		}

		protected override void IncreaseEdgeCapacity(int nEdges)
		{
			base.IncreaseEdgeCapacity(nEdges);

			int oldLength = edgeEvents.Length;
			EdgeEvent[] newEdgeEvents = new EdgeEvent[nEdges];
			Array.Copy(edgeEvents, newEdgeEvents, oldLength);
			edgeEvents = newEdgeEvents;

			for (int i = oldLength; i < maxEdges; ++i)
			{
				edgeEvents[i].eventType = EventType.None;
				edgeEvents[i].queueId = -1;
			}
		}
	}
}