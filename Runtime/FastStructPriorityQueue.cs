using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;

namespace Briganti.StraightSkeletons.Priority_Queue
{
	/// <summary>
	/// An implementation of a min-Priority Queue using a heap.  Has O(1) .Contains()!
	/// See https://github.com/BlueRaja/High-Speed-Priority-Queue-for-C-Sharp/wiki/Getting-Started for more information
	/// </summary>
	/// <typeparam name="T">The values in the queue.  Must extend the FastPriorityQueueNode class</typeparam>
	public sealed class FastStructPriorityQueue<T>
	{
		private struct Node
		{
			public float Priority;
			public int QueueIndex;
			public int Id;

			public override string ToString()
			{
				return $"Node at {QueueIndex} with priority {Priority} and id {Id}";
			}
		}

		private int _maxNodes;

		private int _numNodes;
		private Node[] _nodes;

		private T[] _values;
		private int[] _queueIndexById;
		private List<int> _availableIds;

		/// <summary>
		/// Instantiate a new Priority Queue
		/// </summary>
		/// <param name="maxNodes">The max nodes ever allowed to be enqueued (going over this will cause undefined behavior)</param>
		public FastStructPriorityQueue(int maxNodes)
		{
#if DEBUG
			if (maxNodes <= 0)
			{
				throw new InvalidOperationException("New queue size cannot be smaller than 1");
			}
#endif
			_maxNodes = maxNodes;

			_numNodes = 0;
			_nodes = new Node[maxNodes + 1];

			_values = new T[maxNodes];
			_queueIndexById = new int[maxNodes];

			_availableIds = new List<int>();
			_availableIds.AddRange(Enumerable.Range(0, maxNodes));

			for (int i = 0; i < _queueIndexById.Length; ++i)
			{
				_queueIndexById[i] = -1;
			}
		}

		/// <summary>
		/// Returns the number of nodes in the queue.
		/// O(1)
		/// </summary>
		public int Count
		{
			get
			{
				return _numNodes;
			}
		}

		/// <summary>
		/// Returns the maximum number of items that can be enqueued at once in this queue.  Once you hit this number (ie. once Count == MaxSize),
		/// attempting to enqueue another item will cause undefined behavior.  O(1)
		/// </summary>
		public int MaxSize
		{
			get
			{
				return _nodes.Length - 1;
			}
		}

		/// <summary>
		/// Removes every node from the queue.
		/// O(n) (So, don't do this often!)
		/// </summary>
#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		public void Clear()
		{
			Array.Clear(_nodes, 1, _numNodes);
			_numNodes = 0;

			_availableIds.Clear();
			_availableIds.AddRange(Enumerable.Range(0, _maxNodes));
		}

		/// <summary>
		/// Returns (in O(1)!) whether the given node is in the queue.
		/// If node is or has been previously added to another queue, the result is undefined unless oldQueue.ResetNode(node) has been called
		/// O(1)
		/// </summary>
#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		public bool Contains(int nodeId)
		{
			if (nodeId < 0 || nodeId >= _queueIndexById.Length) return false;
			int queueIndex = _queueIndexById[nodeId];
			return queueIndex != -1;
		}

		/// <summary>
		/// Enqueue a node to the priority queue.  Lower values are placed in front. Ties are broken arbitrarily.
		/// If the queue is full, the result is undefined.
		/// If the node is already enqueued, the result is undefined.
		/// If node is or has been previously added to another queue, the result is undefined unless oldQueue.ResetNode(node) has been called
		/// O(log n)
		/// </summary>
#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		public void Enqueue(T value, float priority, out int nodeId)
		{
#if DEBUG
			if (_numNodes >= _nodes.Length - 1)
			{
				throw new InvalidOperationException("Queue is full - node cannot be added: " + value);
			}
#endif
			nodeId = GetNewNodeId();
			_values[nodeId] = value;

			_numNodes++;
			ref Node node = ref _nodes[_numNodes];
			node.Priority = priority;
			node.QueueIndex = _numNodes;
			node.Id = nodeId;
			_queueIndexById[nodeId] = node.QueueIndex;
			CascadeUp(node);
		}

		private int GetNewNodeId()
		{
			if (_availableIds.Count() == 0) throw new ArgumentException($"Queue is full!");
			int idx = _availableIds.Count - 1;
			int id = _availableIds[idx];
			_availableIds.RemoveAt(idx);
			return id;
		}

#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		private void CascadeUp(Node node)
		{
			//aka Heapify-up
			int parent;
			if (node.QueueIndex > 1)
			{
				parent = node.QueueIndex >> 1;
				ref Node parentNode = ref _nodes[parent];
				if (HasHigherOrEqualPriority(parentNode, node))
					return;

				//Node has lower priority value, so move parent down the heap to make room
				parentNode = UpdateQueueIndex(parentNode, node.QueueIndex);
				node.QueueIndex = parent;
			}
			else
			{
				return;
			}
			while (parent > 1)
			{
				parent >>= 1;
				ref Node parentNode = ref _nodes[parent];
				if (HasHigherOrEqualPriority(parentNode, node))
					break;

				//Node has lower priority value, so move parent down the heap to make room
				parentNode = UpdateQueueIndex(parentNode, node.QueueIndex);
				node.QueueIndex = parent;
			}

			node = UpdateQueueIndex(node, node.QueueIndex);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		private Node UpdateQueueIndex(in Node node, int queueIndex)
		{
			_nodes[queueIndex] = node;
			_queueIndexById[node.Id] = queueIndex;
			_nodes[queueIndex].QueueIndex = queueIndex;
			return _nodes[queueIndex];
		}

#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		private void CascadeDown(Node node)
		{
			//aka Heapify-down
			int finalQueueIndex = node.QueueIndex;
			int childLeftIndex = 2 * finalQueueIndex;

			// If leaf node, we're done
			if (childLeftIndex > _numNodes)
			{
				return;
			}

			// Check if the left-child is higher-priority than the current node
			int childRightIndex = childLeftIndex + 1;
			Node childLeft = _nodes[childLeftIndex];
			if (HasHigherPriority(childLeft, node))
			{
				// Check if there is a right child. If not, swap and finish.
				if (childRightIndex > _numNodes)
				{
					UpdateQueueIndex(node, childLeftIndex);
					UpdateQueueIndex(childLeft, finalQueueIndex);
					return;
				}
				// Check if the left-child is higher-priority than the right-child
				Node childRight = _nodes[childRightIndex];
				if (HasHigherPriority(childLeft, childRight))
				{
					// left is highest, move it up and continue
					UpdateQueueIndex(childLeft, finalQueueIndex);
					finalQueueIndex = childLeftIndex;
				}
				else
				{
					// right is even higher, move it up and continue
					UpdateQueueIndex(childRight, finalQueueIndex);
					finalQueueIndex = childRightIndex;
				}
			}
			// Not swapping with left-child, does right-child exist?
			else if (childRightIndex > _numNodes)
			{
				return;
			}
			else
			{
				// Check if the right-child is higher-priority than the current node
				Node childRight = _nodes[childRightIndex];
				if (HasHigherPriority(childRight, node))
				{
					UpdateQueueIndex(childRight, finalQueueIndex);
					finalQueueIndex = childRightIndex;
				}
				// Neither child is higher-priority than current, so finish and stop.
				else
				{
					return;
				}
			}

			while (true)
			{
				childLeftIndex = 2 * finalQueueIndex;

				// If leaf node, we're done
				if (childLeftIndex > _numNodes)
				{
					UpdateQueueIndex(node, finalQueueIndex);
					break;
				}

				// Check if the left-child is higher-priority than the current node
				childRightIndex = childLeftIndex + 1;
				childLeft = _nodes[childLeftIndex];
				if (HasHigherPriority(childLeft, node))
				{
					// Check if there is a right child. If not, swap and finish.
					if (childRightIndex > _numNodes)
					{
						UpdateQueueIndex(node, childLeftIndex);
						UpdateQueueIndex(childLeft, finalQueueIndex);
						break;
					}
					// Check if the left-child is higher-priority than the right-child
					Node childRight = _nodes[childRightIndex];
					if (HasHigherPriority(childLeft, childRight))
					{
						// left is highest, move it up and continue
						UpdateQueueIndex(childLeft, finalQueueIndex);
						finalQueueIndex = childLeftIndex;
					}
					else
					{
						// right is even higher, move it up and continue
						UpdateQueueIndex(childRight, finalQueueIndex);
						finalQueueIndex = childRightIndex;
					}
				}
				// Not swapping with left-child, does right-child exist?
				else if (childRightIndex > _numNodes)
				{
					UpdateQueueIndex(node, finalQueueIndex);
					break;
				}
				else
				{
					// Check if the right-child is higher-priority than the current node
					Node childRight = _nodes[childRightIndex];
					if (HasHigherPriority(childRight, node))
					{
						UpdateQueueIndex(childRight, finalQueueIndex);
						finalQueueIndex = childRightIndex;
					}
					// Neither child is higher-priority than current, so finish and stop.
					else
					{
						UpdateQueueIndex(node, finalQueueIndex);
						break;
					}
				}
			}
		}

		/// <summary>
		/// Returns true if 'higher' has higher priority than 'lower', false otherwise.
		/// Note that calling HasHigherPriority(node, node) (ie. both arguments the same node) will return false
		/// </summary>
#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		private bool HasHigherPriority(in Node higher, in Node lower)
		{
			return (higher.Priority < lower.Priority);
		}

		/// <summary>
		/// Returns true if 'higher' has higher priority than 'lower', false otherwise.
		/// Note that calling HasHigherOrEqualPriority(node, node) (ie. both arguments the same node) will return true
		/// </summary>
#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		private bool HasHigherOrEqualPriority(in Node higher, in Node lower)
		{
			return (higher.Priority <= lower.Priority);
		}

		/// <summary>
		/// Removes the head of the queue and returns it.
		/// If queue is empty, result is undefined
		/// O(log n)
		/// </summary>
#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		public T Dequeue()
		{
#if DEBUG
			if (_numNodes <= 0)
			{
				throw new InvalidOperationException("Cannot call Dequeue() on an empty queue");
			}
#endif

			Node returnMe = _nodes[1];

			// remember that this id is available again
			RemoveNode(ref returnMe);

			//If the node is already the last node, we can remove it immediately
			if (_numNodes == 1)
			{
				_numNodes = 0;
				return _values[returnMe.Id];
			}

			//Swap the node with the last node
			Node formerLastNode = _nodes[_numNodes];
			formerLastNode = UpdateQueueIndex(formerLastNode, 1);
			_numNodes--;

			//Now bubble formerLastNode (which is no longer the last node) down
			CascadeDown(formerLastNode);
			return _values[returnMe.Id];
		}

		private void RemoveNode(ref Node node)
		{
			_queueIndexById[node.Id] = -1;
			_availableIds.Add(node.Id);
		}

		/// <summary>
		/// Returns the head of the queue, without removing it (use Dequeue() for that).
		/// If the queue is empty, behavior is undefined.
		/// O(1)
		/// </summary>
		public T First
		{
			get
			{
#if DEBUG
				if (_numNodes <= 0)
				{
					throw new InvalidOperationException("Cannot call .First on an empty queue");
				}
#endif

				return _values[_nodes[1].Id];
			}
		}

		/// <summary>
		/// This method must be called on a node every time its priority changes while it is in the queue.  
		/// <b>Forgetting to call this method will result in a corrupted queue!</b>
		/// Calling this method on a node not in the queue results in undefined behavior
		/// O(log n)
		/// </summary>
#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		public void UpdatePriority(int nodeId, float priority)
		{
			var queueIndex = _queueIndexById[nodeId];
			ref Node node = ref _nodes[queueIndex];

			node.Priority = priority;
			OnNodeUpdated(node);
		}

#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		private void OnNodeUpdated(Node node)
		{
			//Bubble the updated node up or down as appropriate
			int parentIndex = node.QueueIndex >> 1;

			if (parentIndex > 0 && HasHigherPriority(node, _nodes[parentIndex]))
			{
				CascadeUp(node);
			}
			else
			{
				//Note that CascadeDown will be called if parentNode == node (that is, node is the root)
				CascadeDown(node);
			}
		}

		/// <summary>
		/// Removes a node from the queue.  The node does not need to be the head of the queue.  
		/// If the node is not in the queue, the result is undefined.  If unsure, check Contains() first
		/// O(log n)
		/// </summary>
#if NET_VERSION_4_5
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
#endif
		public void Remove(int nodeId)
		{
			var queueIndex = _queueIndexById[nodeId];
			ref Node node = ref _nodes[queueIndex];
			RemoveNode(ref node);

			//If the node is already the last node, we can remove it immediately
			if (node.QueueIndex == _numNodes)
			{
				_numNodes--;
				return;
			}

			//Swap the node with the last node
			ref Node formerLastNode = ref _nodes[_numNodes];
			formerLastNode = UpdateQueueIndex(formerLastNode, node.QueueIndex);
			_numNodes--;

			//Now bubble formerLastNode (which is no longer the last node) up or down as appropriate
			OnNodeUpdated(formerLastNode);
		}

		private readonly struct Candidate
		{
			public readonly float Priority;
			public readonly int QueueIndex;

			public Candidate(float priority, int queueIndex)
			{
				Priority = priority;
				QueueIndex = queueIndex;
			}
		}


		public IEnumerable<T> EnumerateFirst(int n)
		{
			if (n < 0 || n > _numNodes) throw new ArgumentOutOfRangeException(nameof(n), $"n must be between 0 and Count (Count={_numNodes}).");
			if (n == 0) yield break;

			// Local min-heap of candidate queue indices (heap nodes), ordered by priority.
			// We never push more than ~2n nodes in typical cases; cap for safety.
			int capacity = Math.Min(_numNodes, Math.Max(4, 2 * n + 1));
			var cand = new Candidate[capacity + 1]; // 1-based heap
			int candCount = 0;

			void EnsureCapacity()
			{
				if (candCount + 1 < cand.Length) return;
				Array.Resize(ref cand, cand.Length * 2);
			}

			void PushCandidate(int queueIndex)
			{
				EnsureCapacity();
				var node = _nodes[queueIndex];

				int i = ++candCount;
				while (i > 1)
				{
					int parent = i >> 1;
					if (cand[parent].Priority <= node.Priority) break;
					cand[i] = cand[parent];
					i = parent;
				}

				cand[i] = new Candidate(node.Priority, queueIndex);
			}

			int PopMinQueueIndex()
			{
				// assumes candCount > 0
				int result = cand[1].QueueIndex;

				Candidate last = cand[candCount--];
				if (candCount == 0) return result;

				int i = 1;
				while (true)
				{
					int left = i << 1;
					if (left > candCount) break;

					int right = left + 1;
					int smallest = (right <= candCount && cand[right].Priority < cand[left].Priority) ? right : left;

					if (cand[smallest].Priority >= last.Priority) break;

					cand[i] = cand[smallest];
					i = smallest;
				}

				cand[i] = last;
				return result;
			}

			// Start from root of your main heap
			PushCandidate(1);

			for (int produced = 0; produced < n; produced++)
			{
				int qi = PopMinQueueIndex();
				yield return _values[_nodes[qi].Id];

				int leftChild = qi << 1;
				if (leftChild <= _numNodes)
				{
					PushCandidate(leftChild);
					int rightChild = leftChild + 1;
					if (rightChild <= _numNodes)
						PushCandidate(rightChild);
				}
			}

			// Local struct for the candidate heap
			// (cannot be declared inside a method in older language versions unless using local functions carefully)
			// So we declare it just below as a private struct in the class (see next snippet).
		}

	}
}