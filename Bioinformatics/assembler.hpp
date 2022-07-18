#pragma once
#include "proccess.h"
#include <boost/config.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>

extern cwd::seqData_t assemblySeq;

namespace cwd
{
	using vertex_iterator = boost::graph_traits<AGraph>::vertex_iterator;
	using edge_iterator = boost::graph_traits<AGraph>::edge_iterator;
	using edge_descriptor = boost::graph_traits<AGraph>::edge_descriptor;
	using vertex_descriptor = boost::graph_traits<AGraph>::vertex_descriptor;
	using out_edge_iterator = boost::graph_traits<AGraph>::out_edge_iterator;
	using in_edge_iterator = boost::graph_traits<AGraph>::in_edge_iterator;
	seqan::Dna5String concatReads(const seqan::Dna5String& pre, uint r2, const seqan::Dna5String& R2, const assemblyInfo_t & ovl, cwd::AEdge::Adj adj, bool shouldRev, bool toPrev);

	void generateContig(const cwd::AGraph & g, std::vector<vertex_descriptor>& p, const cwd::seqData_t& seq, std::ofstream& seqOut, int id, int & totalLen);
	
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

	vector<vector<vertex_descriptor>> findPath(const AGraph & g, const seqData_t& seq, ofstream & outAssembly)
	{
		using namespace boost;
		vertex_iterator vi, vend;
		tie(vi, vend) = vertices(g);
		//uint len = std::distance(vi, vend);
		auto e_prop = boost::get(&AEdge::weight, g);
		auto ovl_prop = boost::get(&AEdge::ovl, g);
		auto adj_prop = boost::get(&AEdge::adj, g);
		auto v_prop = boost::get(vertex_property_t(), g);
		int n = num_vertices(g);
		//vector<pair<vector<int>, size_t>> p_v;
		vector<vector<vertex_descriptor>> p_v;
		outAssembly << " Vn: " << n << endl;
#if 0

		for (size_t src = 0; src < n; src++)
		{
			if (boost::out_degree(src, g) < 1)
			{
				continue;
			}
			vector<double> v_dis(n, numeric_limits<short>::max());
			v_dis[src] = 0;
			vector<int> path(n, -1);
			vector<int> p;

			for (int i = 0; i < n; ++i)
				path[i] = i;
			bool r = boost::bellman_ford_shortest_paths(g, n, e_prop, &path[0], &v_dis[0],
				boost::closed_plus<double>(), std::less<double>(), default_bellman_visitor());
			int tar = std::distance(v_dis.begin(),
				std::min_element(std::next(v_dis.begin()), v_dis.end()));
			if (r)
			{
				//copy(v_dis.begin(), v_dis.end(), ostream_iterator<int, char>(cout, " "));
				//cout << endl;
				//copy(path.begin(), path.end(), ostream_iterator<int, char>(cout, " "));
				//cout << endl;
				outAssembly << " path: ";
				int s;
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
				if (p.size() > 0)
				{
					p_v.push_back({ p, tar });
				}
			}
			else
				cout << "negative cycle." << endl;
		}
		auto p = max_element(p_v.begin(), p_v.end(), 
								[](pair<vector<int>, size_t>& a, pair<vector<int>, size_t>& b)
								{ 
									return a.first.size() < b.first.size();
								}
							)->first;
		vector<vector<int>> ps;
		for (auto& i : p_v)
		{
			if (i.first.size() >= p.size())
			{
				ps.push_back(i.first);
			}
		}
		return ps;
#endif // 0

		out_edge_iterator e, eend;
		std::set<vertex_descriptor> path;
		vector<vertex_descriptor> vex;
		int id = 1;
		//ofstream seqOut("result_ecoli_debug22.fasta", ios_base::app);
		for (auto v = vi; v != vend; v++)
		{
			vex.push_back(*v);
			//auto assembly = seq[v_prop[*v].r];
			//seqOut << boost::format(">contig_%d; length=%d; reads=%d; type=dna\n") % id++ % length(assembly) % 1;
			//auto a = begin(assembly);
			//for (; a < prev(end(assembly), 1000); std::advance(a, 1000))
			//{
			//	seqOut << string{ a, next(a, 1000) } << endl;
			//}
			//seqOut << string{ a, end(assembly) } << endl;
			//seqOut << endl;
		}
		//return p_v;

		for (auto& vs : vex)
		{
			if (path.find(vs) != path.end())
			{
				continue;
			}
			vector<vertex_descriptor> pt;
			vector<vertex_descriptor> pt_pre;
			vertex_descriptor rev_vs = -1;
			vertex_descriptor rev_ve = -1;
			pt.push_back(vs);
			path.insert(vs);
			AEdge::Adj adj;
			bool preOrient = true;
			//cout << v_prop[vs].r << " (start)";
			auto vx = vs;
			auto prev_v = vs;
			auto idx = 0;
			bool flag_dup = false;
			auto cnt = 0;
			while (boost::out_degree(vx, g))
			{
				tie(e, eend) = boost::out_edges(vx, g);
				vector<edge_descriptor> v_edge;
				for (auto ei = e; ei != eend; ei++)
				{
					if (preOrient)
					{
						if (adj == AEdge::HeadTail)
						{
							if (adj_prop[*ei] == AEdge::HeadHead or adj_prop[*ei] == AEdge::HeadTail)
							{
								v_edge.push_back(*ei);
							}
						}
						else if (adj == AEdge::TailHead)
						{
							if (adj_prop[*ei] == AEdge::TailHead or adj_prop[*ei] == AEdge::TailTail)
							{
								v_edge.push_back(*ei);
							}
						}
						else
						{
							v_edge.push_back(*ei);
						}
					}
					else
					{
						if (adj == AEdge::HeadHead)
						{
							if (adj_prop[*ei] == AEdge::TailHead or adj_prop[*ei] == AEdge::TailTail)
							{
								v_edge.push_back(*ei);
							}
						}
						else if (adj == AEdge::TailTail)
						{
							if (adj_prop[*ei] == AEdge::HeadTail or adj_prop[*ei] == AEdge::HeadHead)
							{
								v_edge.push_back(*ei);
							}
						}
						else
						{
							;
						}
					}
				}
				if (v_edge.empty())
				{
					break;
				}
				else
				{
					auto break_flag = false;
					for (auto& e : v_edge)
					{
						auto v = boost::target(e, g);
						if (path.find(v) != path.end())
						{
							cnt++;
							if (find(pt.begin(), pt.end(), v) != pt.end())
							{
								// record v;
								rev_vs = v;
							}
							// record prev;
							rev_ve = vx;
							break;
						}
					}
					if (cnt > 0)
					{
						break_flag = true;
						cnt = 0;
					}
					if (break_flag)
					{
						path.insert(vs);
						break;
					}
				}
				sort(v_edge.begin(), v_edge.end(),
					[&e_prop](edge_descriptor a, edge_descriptor b)
				{
					return e_prop[a] > e_prop[b];
				});
				auto p = v_edge.begin();
				auto orient = ovl_prop[*p].orient;
				auto nextv = boost::target(*p, g);
				while (!path.insert(nextv).second && ((p = std::next(p)) != v_edge.end()))
				{
					flag_dup = true;
					nextv = boost::target(*p, g);
					orient = ovl_prop[*p].orient;
				}
				if (p != v_edge.end())
				{
					auto ovl = ovl_prop[*p];
					auto sp1 = v_prop[vx].r == ovl.r1 ? ovl.SP1 : ovl.SP2;
					auto ep1 = v_prop[vx].r == ovl.r1 ? ovl.EP1 : ovl.EP2;
					string ovl1 = { begin(seq[v_prop[vx].r]) + sp1,  begin(seq[v_prop[vx].r]) + ep1 };
					prev_v = vx;
					vx = nextv;
					pt.push_back(nextv);
					preOrient = orient;
					adj = adj_prop[*p];
					auto sp2 = v_prop[nextv].r == ovl.r2 ? ovl.SP2 : ovl.SP1;
					auto ep2 = v_prop[nextv].r == ovl.r2 ? ovl.EP2 : ovl.EP1;
					string ovl2 = { begin(seq[v_prop[nextv].r]) + sp2,  begin(seq[v_prop[nextv].r]) + ep2 };
					//cout << " -> " << v_prop[nextv].r << " (" << orient
					//	<< boost::format("; adj: %u, SP1: %u; EP1: %u; SP2: %u; EP2: %u; %f")
					//	% adj % sp1 % ep1 % sp2 % ep2 % hamming(ovl1, ovl2)
					//	<< ") " << endl;
				}
				else
				{
					break;
				}
			}
			//cout << endl;
			if (pt.size() > 1)
			{
				vx = pt[0];
				auto f_e = boost::edge(vx, pt[1], g);
				adj = adj_prop[f_e.first];
				preOrient = ovl_prop[f_e.first].orient;
				in_edge_iterator e, eend;
				auto cnt = 0;
				while (boost::in_degree(vx, g))
				{
					tie(e, eend) = boost::in_edges(vx, g);
					vector<edge_descriptor> v_edge;
					for (auto ei = e; ei != eend; ei++)
					{
						if (preOrient)
						{
							if (adj == AEdge::HeadTail)
							{
								if (adj_prop[*ei] == AEdge::TailTail or adj_prop[*ei] == AEdge::HeadTail)
								{
									v_edge.push_back(*ei);
								}
							}
							else if (adj == AEdge::TailHead)
							{
								if (adj_prop[*ei] == AEdge::TailHead or adj_prop[*ei] == AEdge::HeadHead)
								{
									v_edge.push_back(*ei);
								}
							}
							else
							{
								v_edge.push_back(*ei);
							}
						}
						else
						{
							if (adj == AEdge::HeadHead)
							{
								if (adj_prop[*ei] == AEdge::HeadTail or adj_prop[*ei] == AEdge::TailTail)
								{
									v_edge.push_back(*ei);
								}
							}
							else if (adj == AEdge::TailTail)
							{
								if (adj_prop[*ei] == AEdge::TailHead or adj_prop[*ei] == AEdge::HeadHead)
								{
									v_edge.push_back(*ei);
								}
							}
							else
							{
								;
							}
						}
					}
					if (v_edge.empty())
					{
						break;
					}
					else
					{
						auto break_flag = false;
						for (auto& e : v_edge)
						{
							auto v = boost::source(e, g);
							if (path.find(v) != path.end())
							{
								cnt++;
								if (find(pt.begin(), pt.end(), v) != pt.end())
								{
									// record v;
									rev_vs = v;
								}

								else if (find(pt_pre.begin(), pt_pre.end(), v) != pt_pre.end())
								{
									// record v;
									rev_vs = v;
								}

								rev_ve = vx;
								break;
							}
						}
						if (cnt > 0)
						{
							break_flag = true;
							cnt = 0;
						}
						if (break_flag)
						{
							path.insert(vs);
							break;
						}
					}

					sort(v_edge.begin(), v_edge.end(),
						[&e_prop](edge_descriptor a, edge_descriptor b)
					{
						return e_prop[a] > e_prop[b];
					});
					auto p = v_edge.begin();
					auto orient = ovl_prop[*p].orient;
					auto nextv = boost::source(*p, g);
					//cout << vs << " ";
					while (!path.insert(nextv).second && ((p = std::next(p)) != v_edge.end()))
					{
						flag_dup = true;
						nextv = boost::source(*p, g);
						orient = ovl_prop[*p].orient;
					}
					if (p != v_edge.end())
					{
						prev_v = vx;
						vx = nextv;
						pt_pre.push_back(nextv);
						preOrient = orient;
						adj = adj_prop[*p];
						//cout << " -> " << "(" << orient << ") " << nextv;
					}
					else
					{
						break;
					}
				}
				std::reverse(pt_pre.begin(), pt_pre.end());
				pt.insert(pt.begin(), pt_pre.begin(), pt_pre.end());
				auto i_rev_vs = find(pt.begin(), pt.end(), rev_vs);
				auto i_rev_ve = find(pt.begin(), pt.end(), rev_ve);
				if (i_rev_ve != pt.end() and i_rev_vs != pt.end())
				{
					if (i_rev_ve == pt.begin())
					{
						i_rev_ve = pt.end();
						//std::reverse(i_rev_vs, i_rev_ve);
					}
					else
					{
						//std::reverse(std::next(i_rev_vs), std::next(i_rev_ve));
					}
				}
			}

			//cout << endl;
			//in_edge_iterator ei, eiend;
			//vector<vertex_descriptor> pt_pre;
			//while (boost::in_degree(vs, g))
			//{
			//	tie(ei, eiend) = boost::in_edges(vs, g);
			//	auto p = min_element(ei, eiend,
			//		[&e_prop](edge_descriptor a, edge_descriptor b)
			//	{
			//		return e_prop[a] < e_prop[b];
			//	}
			//	);
			//	auto nextv = boost::source(*p, g);
			//	if (source.insert(nextv).second)
			//	{
			//		vs = nextv;
			//		pt_pre.push_back(nextv);
			//		source.insert(nextv);
			//	}
			//	else
			//	{
			//		break;
			//	}
			//}
			//std::reverse(pt_pre.begin(), pt_pre.end());
			//pt.insert(pt.begin(), pt_pre.begin(), pt_pre.end());
			if (pt.size() > 0)
			{
				//std::reverse(pt.begin(), pt.end());
				//copy(pt.begin(), prev(pt.end()), ostream_iterator<int, char>(outAssembly, " -> "));
				//outAssembly << *pt.rbegin() << endl;
				//cout << pt.size() << endl;
				p_v.push_back(pt);
			}

		}
		//auto tmp = p_v[0];
		//p_v.clear();
		//for (auto i = tmp.begin(); i != tmp.end(); i++)
		//{
		//	p_v.push_back({ tmp.begin(), i });
		//}
		//p_v.push_back({ tmp.begin(), tmp.end() });
		return p_v;
		vector<vertex_descriptor> merge_p;
		vector<vector<vertex_descriptor>> p_m;
		for (auto& p : p_v)
		{
			auto v = *p.begin();
			for (auto ip = p_v.begin() + 1; ip != p_v.end(); ip++)
			{
				auto ve = *ip->rbegin();
				auto e = edge(*(p.begin() + 1), v, g);
				auto ei = edge(v, ve, g);
				if (ei.second)
				{
					auto adj = adj_prop[e.first];
					auto adj2 = adj_prop[ei.first];
					if (adj == AEdge::HeadTail)
					{
						if (adj2 == AEdge::HeadHead or adj2 == AEdge::HeadTail)
						{
							merge_p = *ip;
							merge_p.insert(merge_p.end(), p.begin(), p.end());
							p_m.push_back(merge_p);
						}
					}
					else if (adj == AEdge::TailHead)
					{
						if (adj2 == AEdge::TailHead or adj2 == AEdge::TailTail)
						{
							merge_p = *ip;
							merge_p.insert(merge_p.end(), p.begin(), p.end());
							p_m.push_back(merge_p);
						}
					}
					else if (adj == AEdge::HeadHead)
					{
						if (adj2 == AEdge::TailHead or adj2 == AEdge::TailTail)
						{
							merge_p = *ip;
							merge_p.insert(merge_p.end(), p.begin(), p.end());
							p_m.push_back(merge_p);
						}
					}
					else if (adj == AEdge::TailTail)
					{
						if (adj2 == AEdge::HeadTail or adj2 == AEdge::HeadHead)
						{
							merge_p = *ip;
							merge_p.insert(merge_p.end(), p.begin(), p.end());
							p_m.push_back(merge_p);
						}
					}
				}
				else
				{
					p_m.push_back(*ip);
				}
			}
			for (auto ip = p_v.begin() + 1; ip != p_v.end(); ip++)
			{
				auto ve = *ip->begin();
				auto e = edge(*(p.begin() + 1), v, g);
				auto ei = edge(v, ve, g);
				if (ei.second)
				{
					auto adj = adj_prop[e.first];
					auto adj2 = adj_prop[ei.first];
					if (adj == AEdge::HeadTail)
					{
						if (adj2 == AEdge::HeadHead or adj2 == AEdge::HeadTail)
						{
							merge_p = *ip;
							std::reverse(merge_p.begin(), merge_p.end());
							merge_p.insert(merge_p.end(), p.begin(), p.end());
							p_m.push_back(merge_p);
						}
					}
					else if (adj == AEdge::TailHead)
					{
						if (adj2 == AEdge::TailHead or adj2 == AEdge::TailTail)
						{
							merge_p = *ip;
							std::reverse(merge_p.begin(), merge_p.end());
							merge_p.insert(merge_p.end(), p.begin(), p.end());
							p_m.push_back(merge_p);
						}
					}
					else if (adj == AEdge::HeadHead)
					{
						if (adj2 == AEdge::TailHead or adj2 == AEdge::TailTail)
						{
							merge_p = *ip;
							std::reverse(merge_p.begin(), merge_p.end());
							merge_p.insert(merge_p.end(), p.begin(), p.end());
							p_m.push_back(merge_p);
						}
					}
					else if (adj == AEdge::TailTail)
					{
						if (adj2 == AEdge::HeadTail or adj2 == AEdge::HeadHead)
						{
							merge_p = *ip;
							std::reverse(merge_p.begin(), merge_p.end());
							merge_p.insert(merge_p.end(), p.begin(), p.end());
							p_m.push_back(merge_p);
						}
					}
				}
			}
		}

		return p_m;
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

	void generateContig(const cwd::AGraph & g, std::vector<vertex_descriptor>& p, const cwd::seqData_t& seq, std::ofstream& seqOut, int id, int & totalLen)
	{
		bool toPrev = false;
		auto v_prop = boost::get(vertex_property_t(), g);
		ofstream seqOut2("result_dmel_debug_frag.fasta", ios_base::app);

		if (p.size() > 1)
		{
			auto assembly = seq[v_prop[p[0]].r];
			auto E_ovl = get(&AEdge::ovl, g);
			auto E_adj = get(&AEdge::adj, g);
			auto orient = E_ovl[edge(p[0], p[1], g).first].orient;
			bool shouldRev = !orient;

			for (auto i = p.begin(); /*false && */i != prev(p.end()); i++)
			{
				auto e = edge(*i, *std::next(i), g).first;
				auto e2 = edge(*std::next(i), *i, g).first;
				auto src = source(e, g);
				auto des = target(e, g);
				auto b = src == *i && des == *std::next(i);
				auto sr = v_prop[src].r;
				auto dr = v_prop[des].r;
				E_ovl[e].r1;
				E_ovl[e].r2;
				//auto assembly2 = seq[v_prop[*std::next(i)].r];
				//seqOut2 << boost::format(">contig_%d; length=%d; reads=%d; type=dna\n") % id % length(assembly2) % p.size();
				//auto a = begin(assembly2);
				//for (; a < prev(end(assembly2), 1000); std::advance(a, 1000))
				//{
				//	seqOut2 << string{ a, next(a, 1000) } << endl;
				//}
				//seqOut2 << string{ a, end(assembly2) } << endl;
				//seqOut2 << endl;
				//totalLen += length(assembly2);
#if 1
				auto a1 = E_adj[edge(p[0], p[1], g).first];
				auto a2 = E_adj[e2];
				if (a1 == AEdge::HeadTail || a1 == AEdge::TailTail)
				{
					toPrev = true;
				}
				else
				{
					toPrev = false;
				}
				orient = E_ovl[e].orient;
				if (!orient)
				{
					shouldRev = !shouldRev;
					if (i == p.begin())
					{
						assembly = revComp(string{ begin(assembly), end(assembly) });
					}
				}
				if (e.m_eproperty)
				{
					//outAssembly << "PRE:\n" << assembly << " " << endl;
					//outAssembly << "NEWREAD:\n" << seq[v_prop[*i].r] << endl;
					//outAssembly << "OVL:\n" << E_ovl[e].SP1 << " " 
					//	<< E_ovl[e].SP2 << " " << E_ovl[e].EP1 << " " 
					//	<< E_ovl[e].EP2 << endl;
					//outAssembly << "ADJ:\n" << E_adj[e] << endl;
					//cout << "ADJ: " << bitset<2>(E_adj[e]) << " orient: " << E_ovl[e].orient << endl;
					//if (1/*v_prop[*i].r == E_ovl[e].r1*/)
					//{
					//	cout << v_prop[*i].r << "--- SP1: " << E_ovl[e].SP1 << " EP1 : " << E_ovl[e].EP1 << endl
					//		<< v_prop[*std::next(i)].r << "--- SP2: " << E_ovl[e].SP2 << " EP2: " << E_ovl[e].EP2
					//		<< " LEN1: " << length(seq[v_prop[*i].r]) 
					//		<< " LEN2: " << length(seq[v_prop[*std::next(i)].r]) << endl;
					//	
					//	cout << boost::format("%u(%u) -> %u(%u) [r1:%u r2:%u]: orient: %d; adj: %d\n")
					//		% *i % v_prop[*i].r % *(i + 1) % v_prop[*(i + 1)].r % E_ovl[e].r1
					//		% E_ovl[e].r2 % E_ovl[e].orient % E_adj[e];
					//	
					//	compSeqInRange(seq[v_prop[*i].r], seq[v_prop[*std::next(i)].r], v_prop[*i].r,
					//		E_ovl[e].SP1, E_ovl[e].SP2, E_ovl[e].EP1, E_ovl[e].EP2, std::min(E_ovl[e].EP1 - E_ovl[e].SP1,
					//			  E_ovl[e].EP2 - E_ovl[e].SP2), E_ovl[e].orient);
					//}
					//else
					//{
					//	cout << "NEWREAD: SP1: " << E_ovl[e].SP2 << " EP1: " << E_ovl[e].EP2
					//		<< " SP2: " << E_ovl[e].SP1 << " EP2: " << E_ovl[e].EP1
					//		<< " LEN1: " << length(seq[v_prop[*std::next(i)].r])
					//		<< " LEN2: " << length(seq[v_prop[*i].r]) << endl;
					//	
					//	cout << boost::format("%u(%u) -> %u(%u) [r1:%u r2:%u]: orient: %d; adj: %d\n")
					//		% *i % v_prop[*i].r % *(i + 1) % v_prop[*(i + 1)].r % E_ovl[e].r1
					//		% E_ovl[e].r2 % E_ovl[e].orient % E_adj[e];

					//	compSeqInRange(seq[v_prop[*std::next(i)].r], seq[v_prop[*i].r], v_prop[*std::next(i)].r,
					//		E_ovl[e].SP1, E_ovl[e].SP2, E_ovl[e].EP1, E_ovl[e].EP2, std::min(E_ovl[e].EP1 - E_ovl[e].SP1,
					//			E_ovl[e].EP2 - E_ovl[e].SP2), E_ovl[e].orient);
					//}
					
					assembly = concatReads(assembly, v_prop[des].r,
						seq[v_prop[des].r], E_ovl[e], E_adj[e], shouldRev, toPrev);

					//cout << "RES:\n";
					//cout << string{ begin(assembly), begin(assembly) + 60 } + "..." + string{ end(assembly) - 60, end(assembly) } << endl;
					//outAssembly << "RES:\n" << assembly << endl;
				}
				else
				{
					continue;
				}

#endif // 0
			}
			if (true)
			{
				seqOut << boost::format(">contig_%d; length=%d; reads=%d; type=dna\n") % id % length(assembly) % p.size();
				auto a = begin(assembly);
				for (; a < prev(end(assembly), 1000); std::advance(a, 1000))
				{
					seqOut << string{ a, next(a, 1000) } << endl;
				}
				seqOut << string{ a, end(assembly) } << endl;
				seqOut << endl;
				totalLen += length(assembly);
				//appendValue(assemblySeq, assembly);
			}
			//cout << "total: " << length(assembly) << endl;
		}

		else
		{
			auto assembly = seq[v_prop[p[0]].r];
			seqOut << boost::format(">contig_%d; length=%d; reads=%d; type=dna\n") % id % length(assembly) % p.size();
			auto a = begin(assembly);
			for (; a < prev(end(assembly), 1000); std::advance(a, 1000))
			{
				seqOut << string{ a, next(a, 1000) } << endl;
			}
			seqOut << string{ a, end(assembly) } << endl;
			seqOut << endl;
			totalLen += length(assembly);
		}

	}

	seqan::Dna5String concatReads(const seqan::Dna5String& pre, uint r2, const seqan::Dna5String& R2, const assemblyInfo_t & ovl, cwd::AEdge::Adj adj, bool shouldRev, bool toPrev)
	{
		using seqan::length;
		string newRead = "";
		string newRead1 = "";
		string newRead2 = "";
		if (ovl.orient)
		{
			if (adj == AEdge::HeadHead || adj == AEdge::TailTail)
			{
				;
			}
			else if (adj == AEdge::TailHead)
			{
				newRead1 = { begin(pre), end(pre) };
				newRead2 = { begin(R2) + (r2 == ovl.r1 ? ovl.EP1 : ovl.EP2), end(R2) };
				//cout << boost::format("r1: %d, r2: %d, ovl_len: %d\n") % length(newRead1) % length(newRead2) % (ovl.EP1 - ovl.SP1);
			}
			else
			{
				newRead1 = { begin(pre), end(pre) };
				newRead2 = { begin(R2), begin(R2) + (r2 == ovl.r1 ? ovl.SP1 : ovl.SP2)};
				//cout << boost::format("r1: %d, r2: %d, ovl_len: %d\n") % length(newRead1) % length(newRead2) % (ovl.EP1 - ovl.SP1);
			}
		}
		else 
		{
			if (adj == AEdge::HeadHead)
			{
				newRead1 = { begin(pre), end(pre) };
				newRead2 = { begin(R2) + (r2 == ovl.r1 ? ovl.EP1 : ovl.EP2), end(R2) };
				//cout << boost::format("r1: %d, r2: %d, ovl_len: %d\n") % length(newRead1) % length(newRead2) % (ovl.EP1 - ovl.SP1);
			}
			else if (adj == AEdge::TailTail)
			{
				newRead1 = { begin(pre), end(pre) };
				newRead2 = { begin(R2), begin(R2) + (r2 == ovl.r1 ? ovl.SP1 : ovl.SP2) };
				//cout << boost::format("r1: %d, r2: %d, ovl_len: %d\n") % length(newRead1) % length(newRead2) % (ovl.EP1 - ovl.SP1);
			}
			else
			{
				;
			}
		}
		if (shouldRev)
		{
			newRead2 = cwd::revComp(newRead2);
		}
		if (toPrev)
		{
			newRead += newRead2;
			newRead += newRead1;
		}
		else
		{
			newRead += newRead1;
			newRead += newRead2;
		}
		return newRead;
	}
}


//for (auto& pp : p_v)
//{
//	auto t = find(pp.begin(), pp.end(), v);
//	if (t != pp.end())
//	{
//		auto connected = false;
//		auto Len1 = std::distance(t, pp.end());
//		auto Len2 = std::distance(pp.begin(), t);
//		if (Len1 >= pp.size() * 0.2 && Len1 >= pp.size() * 0.2)
//		{
//			break;
//		}
//		if (Len1 > Len2)
//		{
//			// head
//			if (t == pp.end())
//			{
//				break;
//			}
//			auto f_e = boost::edge(*t, *std::next(t), g);
//			if (!f_e.second)
//			{
//				break;
//			}
//			adj = adj_prop[f_e.first];
//			preOrient = ovl_prop[f_e.first].orient;
//			if (preOrient)
//			{
//				if (adj == AEdge::HeadTail)
//				{
//					if (adj_prop[e] == AEdge::TailTail or adj_prop[e] == AEdge::HeadTail)
//					{
//						connected = true;
//					}
//				}
//				else if (adj == AEdge::TailHead)
//				{
//					if (adj_prop[e] == AEdge::TailHead or adj_prop[e] == AEdge::HeadHead)
//					{
//						connected = true;
//					}
//				}
//				else
//				{
//					;
//				}
//			}
//			else
//			{
//				if (adj == AEdge::HeadHead)
//				{
//					if (adj_prop[e] == AEdge::HeadTail or adj_prop[e] == AEdge::TailTail)
//					{
//						connected = true;
//					}
//				}
//				else if (adj == AEdge::TailTail)
//				{
//					if (adj_prop[e] == AEdge::TailHead or adj_prop[e] == AEdge::HeadHead)
//					{
//						connected = true;
//					}
//				}
//				else
//				{
//					;
//				}
//			}
//			if (Len2 > pt.size())
//			{
//				break;
//			}
//			auto pos = pp.erase(pp.begin(), t);
//			std::reverse(pt.begin(), pt.end());
//			pp.insert(pos, pt.begin(), pt.end());
//			break;
//		}
//		else
//		{
//			//tail
//			if (t == pp.begin())
//			{
//				break;
//			}
//			auto f_e = boost::edge(*std::prev(t), *t, g);
//			if (!f_e.second)
//			{
//				break;
//			}
//			adj = adj_prop[f_e.first];
//			preOrient = ovl_prop[f_e.first].orient;
//			if (preOrient)
//			{
//				if (adj == AEdge::HeadTail)
//				{
//					if (adj_prop[e] == AEdge::TailTail or adj_prop[e] == AEdge::HeadTail)
//					{
//						connected = true;
//					}
//				}
//				else if (adj == AEdge::TailHead)
//				{
//					if (adj_prop[e] == AEdge::TailHead or adj_prop[e] == AEdge::HeadHead)
//					{
//						connected = true;
//					}
//				}
//				else
//				{
//					;
//				}
//			}
//			else
//			{
//				if (adj == AEdge::HeadHead)
//				{
//					if (adj_prop[e] == AEdge::HeadTail or adj_prop[e] == AEdge::TailTail)
//					{
//						connected = true;
//					}
//				}
//				else if (adj == AEdge::TailTail)
//				{
//					if (adj_prop[e] == AEdge::TailHead or adj_prop[e] == AEdge::HeadHead)
//					{
//						connected = true;
//					}
//				}
//				else
//				{
//					;
//				}
//			}
//
//			if (Len2 > pt.size())
//			{
//				break;
//			}
//			auto pos = pp.erase(std::next(t), pp.end());
//			pp.insert(pos, pt.begin(), pt.end());
//			break;
//		}
//	}
//}
//break;
