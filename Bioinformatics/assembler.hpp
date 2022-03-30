#pragma once
#include "proccess.h"
#include <boost/config.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/random_spanning_tree.hpp>
namespace cwd
{
	using vertex_iterator = boost::graph_traits<AGraph>::vertex_iterator;
	using edge_iterator = boost::graph_traits<AGraph>::edge_iterator;
	using edge_descriptor = boost::graph_traits<AGraph>::edge_descriptor;
	using vertex_descriptor = boost::graph_traits<AGraph>::vertex_descriptor;
	seqan::Dna5String concatReads(const seqan::Dna5String& pre, uint r2, const seqan::Dna5String& R2, const assemblyInfo_t & ovl, int adj);
	
	void compSeqInRange(const seqan::Dna5String & r1, const seqan::Dna5String & r2, int r, uint start1, uint start2, uint end1, uint end2, uint len, bool orient)
	{
		using namespace std;
		if (len > 30)
		{
			if (!orient)
			{
				string s1 = { begin(r1) + start1, begin(r1) + start1 + len };
				string s2 = cwd::revComp(string{ begin(r2) + end2 - len, begin(r2) + end2 });
				cout << string{ begin(r1),  begin(r1) + 30 } + "...";
				cout << string{ s1.begin(), s1.begin() + 30 } + "..." + string{ s1.end() - 30, s1.end() };
				cout << "..." + string{end(r1) - 30,  end(r1)} << endl;
				cout << string{ begin(r2),  begin(r2) + 30 } + "...";
				cout << string{ s2.begin(), s2.begin() + 30 } + "..." + string{ s2.end() - 30, s2.end() };
				cout << "..." + string{ end(r2) - 30,  end(r2) } << endl;
			}
				
			else
			{
				string s1 = { begin(r1) + start1, begin(r1) + start1 + len };
				string s2 = { begin(r2) + start2, begin(r2) + start2 + len };
				cout << string{ begin(r1),  begin(r1) + 30 } + "...";
				cout << string{ s1.begin(), s1.begin() + 30 } + "..." + string{ s1.end() - 30, s1.end() };
				cout << "..." + string{ end(r1) - 30,  end(r1) } << endl;
				cout << string{ begin(r2),  begin(r2) + 30 } + "...";
				cout << string{ s2.begin(), s2.begin() + 30 } + "..." + string{ s2.end() - 30, s2.end() };
				cout << "..." + string{ end(r2) - 30,  end(r2) } << endl;
			}

		}
	}
	void printPath(vector<int>& distance, vector<int>& path, int n, int start)
	{
		for (auto i = 0; i < n; i++)
		{
			vector<int> p;
			p.push_back(i);
			for (auto s = i; s != start; s = path[s])
			{
				if (s == -1)
				{
					break;
				}
				p.push_back(s);
			}
			if (p.size() == 1)
			{
				p.push_back(i);
			}
			for (auto i = p.rbegin(); i != p.rend(); i++)
			{
				cout << *i << " -> ";
			}
			cout << endl;
		}
	}

	//void bellman_ford(const AGraph& g, int start)
	//{
	//	int n = num_vertices(g);
	//	auto edges = boost::get(AEdge(), g);
	//	edge_iterator e, eend;
	//	tie(e, eend) = boost::edges(g);
	//	vector<int> v_dis(n, 1);
	//	v_dis[start] = 0;
	//	vector<int> path(n, -1);
	//	bool check = false;
	//	for (int i = 0; i < n - 1; i++)
	//	{
	//		//对m条边进行循环
	//		for (auto eit = e; eit != eend; eit++)
	//		{
	//			auto edge = edges[*eit];
	//			auto to = target(*eit, g);
	//			auto from = source(*eit, g);
	//			// 松弛操作
	//			if (v_dis[from] != 1 && v_dis[to] > v_dis[from] + edge.weight)
	//			{
	//				v_dis[to] = v_dis[from] + edge.weight;
	//				path[to] = from;
	//				check = true;
	//			}
	//		}
	//		if (!check)
	//		{
	//			break;  //没有更新，提前退出循环结束算法
	//		}
	//	}
	//	// 3. 检查是否存在负回路（负环）
	//	bool hasCycle = false;
	//	for (auto i = e; i != eend; i++)
	//	{
	//		auto edge = edges[*i];
	//		auto to = target(*i, g);
	//		auto from = source(*i, g);
	//		if (v_dis[from] != 1 && v_dis[to] > v_dis[from] + edge.weight)
	//		{
	//			hasCycle = true;
	//			break;
	//		}
	//	}
	//	if (hasCycle)
	//	{
	//		cout << "存在负权环\n";
	//	}
	//	else
	//	{
	//		auto m = max_element(v_dis.begin(), v_dis.end());
	//		printPath(v_dis, path, n, start);
	//	}
	//}

	void findPath(const AGraph& g, const seqData_t& seq, ofstream & outAssembly)
	{
		using namespace boost;
		vertex_iterator vi, vend;
		tie(vi, vend) = vertices(g);
		uint len = std::distance(vi, vend);
		auto e_prop = boost::get(&AEdge::weight, g);
		auto v_prop = boost::get(vertex_property_t(), g);
		using weightMap = decltype(e_prop);
		int n = num_vertices(g);
		edge_iterator e, eend;
		tie(e, eend) = boost::edges(g);

		vector<int> v_dis(n, numeric_limits<short>::max());
		v_dis[0] = 0;
		vector<int> path(n, -1);
		for (int i = 0; i < n; ++i)
			path[i] = i;
		bool r = boost::bellman_ford_shortest_paths(g, n, e_prop, &path[0], &v_dis[0],
			boost::closed_plus<int>(), std::less<int>(), default_bellman_visitor());
		int tar = std::distance(v_dis.begin(), 
			std::min_element(std::next(v_dis.begin()), v_dis.end()));
		ofstream seqOut("result.txt", ios_base::app);
		if (r)
		{
			//copy(v_dis.begin(), v_dis.end(), ostream_iterator<int, char>(cout, " "));
			//cout << endl;
			//copy(path.begin(), path.end(), ostream_iterator<int, char>(cout, " "));
			//cout << endl;
			outAssembly << " Vn: " << n << ", path: ";
			int s;
			vector<int> p;
			p.push_back(tar);
			for (s = path[tar]; s != path[s] && s != tar; s = path[s])
			{
				p.push_back(s);
			}
			if (s != tar)
			{
				p.push_back(s);
			}
			std::reverse(p.begin(), p.end());
			copy(p.begin(), prev(p.end()), ostream_iterator<int, char>(outAssembly, " -> "));
			outAssembly << *p.rbegin() << endl;
			if (p.size() > 1)
			{
				auto assembly = seq[v_prop[p[0]].r];
				auto E_ovl = get(&AEdge::ovl, g);
				auto E_adj = get(&AEdge::adj, g);
				for (auto i = p.begin(); /*false && */i != prev(p.end()); i++)
				{
					auto e = edge(*i, *std::next(i), g).first;
					if (e.m_eproperty)
					{
						//outAssembly << "PRE:\n" << assembly << " " << endl;
						//outAssembly << "NEWREAD:\n" << seq[v_prop[*i].r] << endl;
						//outAssembly << "OVL:\n" << E_ovl[e].SP1 << " " 
						//	<< E_ovl[e].SP2 << " " << E_ovl[e].EP1 << " " 
						//	<< E_ovl[e].EP2 << endl;
						//outAssembly << "ADJ:\n" << E_adj[e] << endl;
						//cout << "ADJ: " << bitset<2>(E_adj[e]) << " orient: " << E_ovl[e].orient << endl;
						if (v_prop[*i].r == E_ovl[e].r1)
						{
							//cout << "NEWREAD: SP1: " << E_ovl[e].SP1 << " EP1: " << E_ovl[e].EP1 
							//	<< " SP2: " << E_ovl[e].SP2 << " EP2: " << E_ovl[e].EP2 
							//	<< " LEN1: " << length(seq[v_prop[*i].r]) 
							//	<< " LEN2: " << length(seq[v_prop[*std::next(i)].r]) << endl;
							//compSeqInRange(seq[v_prop[*i].r], seq[v_prop[*std::next(i)].r], v_prop[*i].r,
							//	E_ovl[e].SP1, E_ovl[e].SP2, E_ovl[e].EP1, E_ovl[e].EP2, std::min(E_ovl[e].EP1 - E_ovl[e].SP1,
							//		  E_ovl[e].EP2 - E_ovl[e].SP2), E_ovl[e].orient);
						}
						else
						{
							//cout << "NEWREAD: SP1: " << E_ovl[e].SP2 << " EP1: " << E_ovl[e].EP2
							//	<< " SP2: " << E_ovl[e].SP1 << " EP2: " << E_ovl[e].EP1
							//	<< " LEN1: " << length(seq[v_prop[*std::next(i)].r])
							//	<< " LEN2: " << length(seq[v_prop[*i].r]) << endl;

							//compSeqInRange(seq[v_prop[*std::next(i)].r], seq[v_prop[*i].r], v_prop[*std::next(i)].r,
							//	E_ovl[e].SP1, E_ovl[e].SP2, E_ovl[e].EP1, E_ovl[e].EP2, std::min(E_ovl[e].EP1 - E_ovl[e].SP1,
							//		E_ovl[e].EP2 - E_ovl[e].SP2), E_ovl[e].orient);
						}
						assembly = concatReads(assembly, v_prop[*std::next(i)].r,
							seq[v_prop[*std::next(i)].r], E_ovl[e], E_adj[e]);
						//cout << "RES:\n";
						//cout << string{ begin(assembly), begin(assembly) + 60 } + "..." + string{ end(assembly) - 60, end(assembly) } << endl;
						//outAssembly << "RES:\n" << assembly << endl;
					}
					else
					{
						continue;
						//auto e = edge(*i, *std::next(i), g).first;
						//if (e.m_eproperty)
						//{
						//	cout << "ADJ: " << E_adj[e] << endl;
						//	cout << "NEWREAD:\n";
						//	compSeqInRange(seq[v_prop[*i].r], seq[v_prop[*std::next(i)].r], v_prop[*std::next(i)].r,
						//		E_ovl[e].SP1, E_ovl[e].SP2, E_ovl[e].EP1, E_ovl[e].EP2, std::min(E_ovl[e].EP1 - E_ovl[e].SP1,
						//			E_ovl[e].EP2 - E_ovl[e].SP2), E_ovl[e].orient);
						//	assembly = concatReads(assembly, v_prop[*std::next(i)].r,
						//		seq[v_prop[*std::next(i)].r], E_ovl[e], E_adj[e]);
						//	cout << "RES:\n";
						//	cout << string{ begin(assembly), begin(assembly) + 60 } + "..." + string{ end(assembly) - 60, end(assembly) } << endl;
					}
					seqOut << length(assembly) << endl;
					auto a = begin(assembly);
					for (auto b = next(a, 1000); b < end(assembly); std::advance(a, 1000), std::advance(b, 1000))
					{
						seqOut << string{ a, b } << endl;
					}
					seqOut << endl;
				}
				//cout << "total: " << length(assembly) << endl;
			}
		}
		else
			cout << "negative cycle." << endl;

		//vector<vertex_descriptor> path;
		
		//using aedge =  property_traits<weightMap>::value_type;
		//int d[450][450];
		//int paths[450][450];
		////auto compare = [](AEdge a, AEdge b) { return a.adj % 2 == b.adj % 2;};
		////auto combine = [](AEdge a, AEdge b) { return a.adj + b.adj;};
		////boost::floyd_warshall_all_pairs_shortest_paths(g, d, e_prop, compare, combine, std::numeric_limits<int>::max(), 0);
		//for (int i = 0; i < 450; i++)
		//{
		//	for (int j = 0; j < 450; j++)
		//	{
		//		paths[i][j] = j;
		//		d[i][j] = 0;
		//	}
		//}
		//for (vertex_descriptor i = *vi; i < *vend; i++)
		//{
		//	for (vertex_descriptor j = *std::next(vi); j < *vend; j++)
		//	{
		//		auto e = edge(i, j, g).first;
		//		if (e.m_eproperty)
		//		{
		//			d[i][j] = e_prop[e].weight;
		//			d[j][i] = d[i][j];
		//		}
		//		else
		//		{
		//			d[i][j] = d[j][i] = 1;
		//		}
		//	}
		//}
		//for (vertex_descriptor k = *vi; k < *vend; k++)
		//{
		//	for (vertex_descriptor i = *vi; i < *vend; i++)
		//	{
		//		for (vertex_descriptor j = *vi; j < *vend; j++)
		//		{
		//			if (d[i][j] > d[i][k] + d[k][j])
		//			{
		//				auto e = edge(i, k, g).first;
		//				auto e2 = edge(k, j, g).first;
		//				if (e.m_eproperty && e2.m_eproperty && e_prop[e].adj % 2 != e_prop[e2].adj % 2)
		//				{
		//					d[i][j] = d[i][k] + d[k][j];
		//					paths[i][j] = k;
		//				}
		//			}
		//		}
		//	}
		//}
		//
		//int max = 0;
		//int scr = 0, des = 0;
		//for (int i = 0; i < len; i++)
		//{
		//	for (int j = 0; j < len; j++)
		//	{
		//		if (max > d[i][j])
		//		{
		//			max = d[i][j];
		//			des = j;
		//			scr = i;
		//		}
		//	}
		//}
		//
		//int next = paths[scr][des];
		//cout << "scr: " << scr << " des: " << des << endl;
		//cout << "path: " << scr << "--" << next << "--";
		//int i = 0;
		//for (i = 0; i < len; i++)
		//{
		//	next = paths[next][des];
		//	if (next == des)
		//	{
		//		break;
		//	}
		//	cout << next << "--";
		//}
		//if (i != len)
		//{
		//	cout << des << endl;
		//}
		
		//for (vertex_descriptor n = *vi; n != *vend; )
		//{
		//	tie(eit, eend) = out_edges(n, g);
		//	if (path.empty())
		//		path.push_back(n);
		//	s[n] = true;
		//	auto pre = e_prop[*eit];
		//	int c = std::distance(eit, eend);
		//	if (c == 0)
		//	{
		//		break;
		//	}
		//	edge_iterator i = eit;

		//	for (i = eit; i != eend; i++)
		//	{
		//		//std::cout << boost::target(it, g) << '\n'; 
		//		edge_descriptor it = *i;
		//		if (e_prop[it].adj % 2 != pre.adj % 2)
		//		{
		//			s[n] = true;
		//			pre.adj = e_prop[it].adj;
		//			if (s[target(it, g)] == true)
		//			{
		//				goto L;
		//			}
		//			else
		//			{
		//				n = target(it, g);
		//				path.push_back(n);
		//				break;
		//			}
		//		}
		//		else
		//		{
		//			;
		//		}
		//	}
		//	if (i == eend)
		//	{
		//		break;
		//	}
		//}




	}

	seqan::Dna5String concatReads(const seqan::Dna5String& pre, uint r2, const seqan::Dna5String& R2, const assemblyInfo_t & ovl, int adj)
	{
		string newRead = "";
		if (ovl.orient)
		{
			if (adj == AEdge::HeadHead || adj == AEdge::TailTail)
			{
				return length(pre) > length(R2) ? pre : R2;
			}
			else if (adj == AEdge::TailHead)
			{
				string newRead1 = { begin(pre), end(pre) };
				string newRead2 = { begin(R2) + (r2 == ovl.r1 ? ovl.EP1 : ovl.EP2), end(R2) };
				newRead += newRead1;
				newRead += newRead2;
			}
			else
			{
				string newRead1 = { begin(R2), begin(R2) + (r2 == ovl.r1 ? ovl.SP1 : ovl.SP2) };
				string newRead2 = { begin(pre), end(pre) };
				newRead += newRead1;
				newRead += newRead2;
			}
		}
		else 
		{
			if (adj == AEdge::HeadHead || adj == AEdge::TailTail)
			{
				return length(pre) > length(R2) ? pre : R2;
			}
			else if (adj == AEdge::TailHead)
			{
				string newRead1 = { begin(pre), end(pre) };
				string newRead2 = { begin(R2), begin(R2) + (r2 == ovl.r1 ? ovl.SP1 : ovl.SP2) };
				newRead += newRead1;
				newRead += cwd::revComp(newRead2);
			}
			else
			{
				string newRead1 = { begin(R2) + (r2 == ovl.r1 ? ovl.SP1 : ovl.SP2) , end(R2)};
				string newRead2 = { begin(pre), end(pre) };
				newRead += cwd::revComp(newRead1);
				newRead += newRead2;
			}
		}
		return newRead;
	}
}


