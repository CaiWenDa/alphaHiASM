#include "assembly.h"
#include "utility.h"
#include "boost/graph/subgraph.hpp"
#include "boost/graph/graph_utility.hpp"
#include "boost/graph/connected_components.hpp"
#include "boost/property_map/property_map.hpp"
#include <iostream>
#include <algorithm>
#include <random>
#include <cstdlib>

#define DEBUG_ONE_BY_ONE 0

using namespace cwd;
using namespace std;

extern cwd::seqData_t assemblySeq;
extern vector<cwd::assemblyInfo_t> overlap;

vector<cwd::ComponentGraph> comps;
shared_ptr<cwd::AGraph> assemblyGraph;
shared_ptr<cwd::SubGraph> tmpGraph;
int OVL_TIP_LEN = 200;// * 80;
int genomeSize = 0;//100 * 1000000 * 1.0; //5000000;
std::set<size_t> delReads;
vector<shared_ptr<list<cwd::assemblyInfo_t>>> assemblyChain;


uint cwd::heuristic(edge_descriptor & p, const cwd::AGraph& g, boost::adj_list_edge_property_map<boost::bidirectional_tag, cwd::AEdge::Adj, const cwd::AEdge::Adj&, cwd::vertex_descriptor, const cwd::AEdge, cwd::AEdge::Adj cwd::AEdge::*>& adj_prop, boost::adj_list_edge_property_map<boost::bidirectional_tag, cwd::assemblyInfo_t, const cwd::assemblyInfo_t&, cwd::vertex_descriptor, const cwd::AEdge, cwd::assemblyInfo_t cwd::AEdge::*>& ovl_prop)
{
	auto nextv = boost::source(p, g);
	auto out_deg = boost::out_degree(nextv, g);
	auto adj = adj_prop[p];
	auto tip_len = ovl_prop[p].SP1;
	if (adj == AEdge::HeadHead or adj == AEdge::HeadTail)
	{
		tip_len = ovl_prop[p].SP1;
	}
	else
	{
		tip_len = ovl_prop[p].SP2;
	}
	uint weight = out_deg * 10 - tip_len;
	return weight * 10;
}

vector<vector<vertex_descriptor>> cwd::findPath(const AGraph & g, const seqData_t& seq, ofstream & outAssembly)
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
		size_t path_id = 0;
		auto pp = p_v.end();
		auto cnt = 0;
		auto pre_ins = true;
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
				auto vv = 0lu;
				for (auto& e : v_edge)
				{
					auto v = boost::target(e, g);
					if (path.find(v) != path.end())
					{
						cnt++;
						vv = v;
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
					for (auto ip = p_v.begin(); ip != p_v.end(); ++ip)
					{
						if (find(ip->begin(), ip->end(), vv) != ip->end())
						{
							pp = ip;
							break;
						}
					}
					break;
				}
			}
			//TODO: A-Star:
			sort(v_edge.begin(), v_edge.end(),
				[&e_prop](edge_descriptor a, edge_descriptor b)
			{
				return e_prop[a] > e_prop[b];
			});
#if DEBUG_ONE_BY_ONE
			ofstream seqOut2("other_reads.fasta");
			for (auto& e : v_edge)
			{
				auto v = boost::target(e, g);
				auto assembly2 = seq[v_prop[v].r];
				seqOut2 << boost::format(">contig_%d; length=?; reads=?; type=dna\n") % ii++;
				auto a = begin(assembly2);
				for (; a < prev(end(assembly2), 1000); std::advance(a, 1000))
				{
					seqOut2 << string{ a, next(a, 1000) } << endl;
				}
				seqOut2 << string{ a, end(assembly2) } << endl;
				seqOut2 << endl;
			}
#endif
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
				//auto ovl = ovl_prop[*p];
				//auto sp1 = v_prop[vx].r == ovl.r1 ? ovl.SP1 : ovl.SP2;
				//auto ep1 = v_prop[vx].r == ovl.r1 ? ovl.EP1 : ovl.EP2;
				//string ovl1 = { begin(seq[v_prop[vx].r]) + sp1,  begin(seq[v_prop[vx].r]) + ep1 };
				prev_v = vx;
				vx = nextv;
				pt.push_back(nextv);
				preOrient = orient;
				adj = adj_prop[*p];
				//auto sp2 = v_prop[nextv].r == ovl.r2 ? ovl.SP2 : ovl.SP1;
				//auto ep2 = v_prop[nextv].r == ovl.r2 ? ovl.EP2 : ovl.EP1;
				//string ovl2 = { begin(seq[v_prop[nextv].r]) + sp2,  begin(seq[v_prop[nextv].r]) + ep2 };
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
			vector<vertex_descriptor> pt_pre;
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
					auto vv = 0lu;
					for (auto& e : v_edge)
					{
						auto v = boost::source(e, g);
						if (path.find(v) != path.end())
						{
							cnt++;
							vv = v;
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
						for (auto ip = p_v.begin(); ip != p_v.end(); ++ip)
						{
							if (find(ip->begin(), ip->end(), vv) != ip->end())
							{
								pp = ip;
								break;
							}
						}
						pre_ins = false;
						break;
					}
				}

				//sort(v_edge.begin(), v_edge.end(),
				//	[&e_prop](edge_descriptor a, edge_descriptor b)
				//{
				//	return e_prop[a] > e_prop[b];
				//});
				auto p = max_element(v_edge.begin(), v_edge.end(),
					[&e_prop, &g, &adj_prop, &ovl_prop](edge_descriptor a, edge_descriptor b)
				{
					return e_prop[a] + heuristic(a, g, adj_prop, ovl_prop) < e_prop[b] + heuristic(b, g, adj_prop, ovl_prop);
				});
#if DEBUG_ONE_BY_ONE
				ofstream seqOut2("other_reads.fasta", ios_base::app);
				for (auto& e : v_edge)
				{
					auto v = boost::target(e, g);
					auto assembly2 = seq[v_prop[v].r];
					seqOut2 << boost::format(">contig_%d; length=?; reads=?; type=dna\n") % ii++;
					auto a = begin(assembly2);
					for (; a < prev(end(assembly2), 1000); std::advance(a, 1000))
					{
						seqOut2 << string{ a, next(a, 1000) } << endl;
					}
					seqOut2 << string{ a, end(assembly2) } << endl;
					seqOut2 << endl;
				}
#endif

				//auto p = v_edge.begin();
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
			++path_id;
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
		if (pt.size() > 1)
		{
			//std::reverse(pt.begin(), pt.end());
			//copy(pt.begin(), prev(pt.end()), ostream_iterator<int, char>(outAssembly, " -> "));
			//outAssembly << *pt.rbegin() << endl;
			//cout << pt.size() << endl;
			if (pp != p_v.end() && pp->size() > 2)
			{
				if (pre_ins)
				{
					//p_v.insert(pp, pt);
					//pp->insert(pp->end(), -1);
					pp->insert(pp->end(), pt.begin(), pt.end());
				}
				else
				{
					//p_v.insert(std::next(pp), pt);
					//pp->insert(pp->begin(), -1);
					pp->insert(pp->begin(), pt.begin(), pt.end());
				}
			}
			else
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

void cwd::generateContig(const cwd::AGraph & g, std::vector<vertex_descriptor>& p, const cwd::seqData_t& seq, std::ofstream& seqOut, int id, int & totalLen)
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
		auto a1 = E_adj[edge(p[0], p[1], g).first];

		for (auto i = p.begin(); /*false && */i != prev(p.end()); i++)
		{
			auto es = edge(*i, *std::next(i), g);
			auto es2 = edge(*std::next(i), *i, g);
			auto e = es.first;
			auto e2 = es2.first;
			auto src = source(e, g);
			auto des = target(e, g);
			auto b = src == *i && des == *std::next(i);
			auto sr = v_prop[src].r;
			auto dr = v_prop[des].r;
			if (es.second && es2.second)
			{
#if 1
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
			else
			{
				//advance(i, 2);
				assembly = concatReadsDirect(assembly, seq[v_prop[*next(i)].r]);
				a1 = E_adj[edge(*next(i), *next(i, 2), g).first];
				shouldRev = !E_ovl[edge(*next(i), *next(i, 2), g).first].orient;
			}
			//E_ovl[e].r1;
			//E_ovl[e].r2;
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

seqan::Dna5String cwd::concatReads(const seqan::Dna5String& pre, uint r2, const seqan::Dna5String& R2, const assemblyInfo_t & ovl, cwd::AEdge::Adj adj, bool shouldRev, bool toPrev)
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

seqan::Dna5String cwd::concatReadsDirect(const seqan::Dna5String& pre, const seqan::Dna5String& R2)
{
	string pr = { begin(pre), end(pre) };
	string r2 = { begin(R2), end(R2) };
	return pr + r2;
}

void cwd::connected_components_subgraphs(AGraph const& g)
{
	using namespace boost;
	boost::copy_graph(g, *tmpGraph);
	std::shared_ptr<vector<unsigned long>> mapping = std::make_shared<vector<unsigned long>>(num_vertices(g));
	size_t num = connected_components(*tmpGraph, mapping->data());

	auto& component_graphs = comps;

	for (size_t i = 0; i < num; i++)
		component_graphs.emplace_back(g,
			[mapping, i, &g](AGraph::edge_descriptor e) {
		return mapping->at(source(e, g)) == i
			|| mapping->at(target(e, g)) == i;
	},
			[mapping, i](AGraph::vertex_descriptor v) {
		return mapping->at(v) == i;
	});

	//return component_graphs;
}

void cwd::createOverlapGraph(seqData_t& seq, size_t block1, size_t block2)
{
	using vertex_descriptor = boost::graph_traits<AGraph>::vertex_descriptor;
	vertex_descriptor src, dst;
	uint idx = 0;
	assemblyGraph = make_shared<AGraph>();
	tmpGraph = make_shared<SubGraph>();
	auto e_prop = boost::get(&AEdge::weight, *assemblyGraph);
	for (auto& ovl : overlap)
	{
		//auto& ovl = overlap[i];
		if (true)
		{
			//int j = length(seq[i]);
			bool flag = false;
			uint r = ovl.r1;
			uint i = ovl.r2;
			uint len_r = length(seq[r]);
			uint len_i = length(seq[i]);
			bool isHeadTail = false, isTailHead = false, isHeadHead = false, isTailTail = false;
			if (
				(isTailHead = len_r - ovl.EP1 < OVL_TIP_LEN && ovl.SP2 < OVL_TIP_LEN && ovl.orient) ||
				(isHeadTail = len_i - ovl.EP2 < OVL_TIP_LEN && ovl.SP1 < OVL_TIP_LEN && ovl.orient) ||
				(isHeadHead = ovl.SP1 < OVL_TIP_LEN && ovl.SP2 < OVL_TIP_LEN && !ovl.orient) ||
				(isTailTail = len_r - ovl.EP1 < OVL_TIP_LEN && len_i - ovl.EP2 < OVL_TIP_LEN && !ovl.orient)
				) //TODO direction
			{
				//delReads.insert(r);
				//delReads.insert(i);
				if (0)
				{
					for (auto& chain : assemblyChain)
					{
						auto ix = find_if(chain->begin(), chain->end(),
							[=](assemblyInfo_t& a) {
							return a.r1 == r || a.r1 == i || a.r2 == r || a.r2 == i;
						});
						if (ix != chain->end())
						{
							//chain->push_back(ovl);
							flag = true;
							chain->insert(next(ix), ovl);
							//break;
						}
						//if (chain->rbegin()->r1 == r || chain->rbegin()->r1 == i)
						//{
						//	chain->push_back({ r, i, ovl.SP1, ovl.EP1, ovl.SP2, ovl.EP2, true });
						//	flag = true;
						//	break;
						//}
						//else if (chain->rbegin()->r2 == r || chain->rbegin()->r2 == i)
						//{
						//	chain->push_back({ i, r, ovl.SP1, ovl.EP1, ovl.SP2, ovl.EP2, true });
						//	flag = true;
						//	break;
						//}
					}
					if (!flag)
					{
						auto n_chain = make_shared<list<assemblyInfo_t>>();
						n_chain->push_back(ovl);
						assemblyChain.push_back(n_chain);
					}
				}
				string ovl1 = { begin(seq[r]) + ovl.SP1,  begin(seq[r]) + ovl.EP1 };
				string ovl2 = { begin(seq[i]) + ovl.SP2,  begin(seq[i]) + ovl.EP2 };
				if (!ovl.orient)
				{
					ovl2 = cwd::revComp(ovl2);
				}
				auto precision = 1.0f - hamming(ovl1, ovl2);
				if (precision < 0.99f)
				{
					continue;
				}
				//if (delReads.find(r) != delReads.end() || delReads.find(i) != delReads.end())
				//{
				//	continue;
				//}
				double weight = 1.0 * (len_r + len_i - ovl.EP1 + ovl.SP1) * precision;
				//double weight = 1.0 * (ovl.EP1 - ovl.SP1) * precision;

				auto& aGraph = assemblyGraph;
				boost::graph_traits<AGraph>::vertex_iterator vi, vi_end;
				boost::tie(vi, vi_end) = boost::vertices(*aGraph);
				auto prop = boost::get(vertex_property_t(), *aGraph);
				auto vx = find_if(vi, vi_end,
					[=](vertex_descriptor ix) {
					AVertex y = prop[ix];
					return y.r == r || y.r == i;
				});
				if (vx != vi_end)
				{
					flag = true;
					AVertex v = { ovl.r1, ovl.SP1, ovl.EP1 };
					AVertex v2 = { ovl.r2, ovl.SP2, ovl.EP2 };
					if (prop[*vx].r == r) // if vx == r1
					{
						auto vx2 = find_if(vi, vi_end,
							[=](vertex_descriptor ix) {
							AVertex y = prop[ix];
							return y.r == i;
						});
						if (vx2 != vi_end)
						{
							auto e = boost::edge(*vx, *vx2, *aGraph);
							if (e.second)
							{
								auto edg = e.first;
								e_prop[edg] += weight;
								//TODO
								auto d = target(edg, *aGraph);
								auto s = source(edg, *aGraph);
								auto p = boost::out_edges(s, *aGraph);
								for (auto i = p.first; i != p.second; i++)
								{
									auto x = target(*i, *aGraph);
									if (x == d)
									{
										e_prop[edg] += weight;
										break;
									}
								}
								if (0/*目的节点在源点的其他边上*/)
								{
									/*e_prop[edg] = weight;*/
								}
							}

							else
							{
								if (isTailHead)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
								}
								else if (isHeadTail)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								}
								else if (isHeadHead)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
								}
								else if (isTailTail)
								{
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
								}
								else
								{
									cerr << "none\n";
								}

							}
							//add to an overlap vector
						}
						else
						{
							src = boost::add_vertex(v2, *aGraph);
							// src == r2
							if (isHeadTail)
							{
								boost::add_edge(*vx, src, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
								// *vx -> src
								boost::add_edge(src, *vx, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								// src -> *vx
							}
							else if (isHeadHead)
							{
								boost::add_edge(*vx, src, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
								boost::add_edge(src, *vx, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
							}
							else if (isTailHead)
							{
								boost::add_edge(*vx, src, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								boost::add_edge(src, *vx, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
							}
							else if (isTailTail)
							{
								boost::add_edge(*vx, src, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
								boost::add_edge(src, *vx, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
							}
							else
							{
								cerr << "none\n";
							}
							//TODO: condition of direction
						}
					}
					else // vx == r2
					{
						auto vx2 = find_if(vi, vi_end,
							[=](vertex_descriptor ix) {
							AVertex y = prop[ix];
							return y.r == r;
						});
						if (vx2 != vi_end)
						{
							auto e = boost::edge(*vx, *vx2, *aGraph);
							if (e.first.m_eproperty)
							{
								auto edg = e.first;
								e_prop[edg] += weight;
								//if (e_prop[edg] < weight)
								//{
									//e_prop[edg] = weight;
								//}
							}

							else
							{
								if (isTailHead)
								{
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
								}
								else if (isHeadTail)
								{
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								}
								else if (isHeadHead)
								{
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
								}
								else if (isTailTail)
								{
									boost::add_edge(*vx2, *vx, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
									boost::add_edge(*vx, *vx2, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
								}
								else
								{
									cerr << "none\n";
								}
							}
							//add to an overlap vector
						}
						else
						{
							dst = boost::add_vertex(v, *aGraph);
							// dst == r1
							if (isTailHead)
							{
								boost::add_edge(dst, *vx, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
								// dst -> *vx
								boost::add_edge(*vx, dst, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
								// *vx -> dst
							}
							else if (isHeadTail)
							{
								boost::add_edge(dst, *vx, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
								boost::add_edge(*vx, dst, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
							}
							else if (isHeadHead)
							{
								boost::add_edge(dst, *vx, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
								boost::add_edge(*vx, dst, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
							}
							else if (isTailTail)
							{
								boost::add_edge(dst, *vx, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
								boost::add_edge(*vx, dst, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
							}
							else
							{
								cerr << "none\n";
							}
							//TODO: condition of direction
						}
					}
				}

				if (!flag)
				{
					auto& aGraph = assemblyGraph;
					AVertex v = { ovl.r1, ovl.SP1, ovl.EP1 };
					AVertex v2 = { ovl.r2, ovl.SP2, ovl.EP2 };
					src = boost::add_vertex(v, *aGraph);
					dst = boost::add_vertex(v2, *aGraph);
					if (isTailHead)
					{
						boost::add_edge(src, dst, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
						boost::add_edge(dst, src, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
					}
					else if (isHeadTail)
					{
						boost::add_edge(src, dst, AEdge{ AEdge::HeadTail, ovl, weight }, *aGraph);
						boost::add_edge(dst, src, AEdge{ AEdge::TailHead, ovl, weight }, *aGraph);
					}
					else if (isHeadHead)
					{
						boost::add_edge(src, dst, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
						boost::add_edge(dst, src, AEdge{ AEdge::HeadHead, ovl, weight }, *aGraph);
					}
					else if (isTailTail)
					{
						boost::add_edge(src, dst, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
						boost::add_edge(dst, src, AEdge{ AEdge::TailTail, ovl, weight }, *aGraph);
					}
					else
					{
						cerr << "none\n";
					}
				}
			}
		}
	}
	connected_components_subgraphs(*assemblyGraph);
}

void cwd::assembler(const seqData_t& seq, ofstream& seqOut)
{
	//ofstream outAssembly("toyAssembly.csv");
	//for (auto& chain : assemblyChain)
	//{
	//	if (chain->size() < 3)
	//	{
	//		continue;
	//	}
	//	for (auto& ovl : *chain)
	//	{
	//		outAssembly << boost::format("%u, %u, %u, %u, %u, %u\n") 
	//			% ovl.r1 % ovl.r2 % ovl.SP1 % ovl.EP1 % ovl.SP2 % ovl.EP2;
	//		
	//	}
	//	//cout << "--------------------------------------------------------------\n";
	//	outAssembly << endl;
	//}
	//assemblyChain.clear();
	//outAssembly.close();

	ofstream outAssembly, outPath;
	outPath.open("toyAssembly_path.txt");
	outAssembly.open("toyAssembly_graph.txt");
	int i = 0;
	int totalLen = 0;
	for (auto& aGraph : comps)
	{
		//	boost::dynamic_properties dp;
		//	dp.property("node_id", boost::get(&AVertex::r, *aGraph));
		//	dp.property("node_id", boost::get(&AVertex::r, *aGraph));
		//	dp.property("label", boost::get(vertex_property_t(), *aGraph));
		//	dp.property("label", boost::get(edge_property_t(), *aGraph));
		//	boost::write_graphviz(outAssembly, *aGraph, dp);
		AGraph g;
		copy_graph(aGraph, g);
		boost::write_graphviz(outAssembly, g, boost::make_label_writer(boost::get(vertex_property_t(), aGraph)),
			make_edge_writer(boost::get(&AEdge::adj, g), boost::get(&AEdge::weight, g)));
		outAssembly << endl;
		outPath << "Graph: " << i++;
		auto ps = findPath(g, seq, outPath);
		sort(ps.begin(), ps.end(),
			[](vector<vertex_descriptor>& a, vector<vertex_descriptor>& b)
		{
			return a.size() > b.size();
		});
		float sum = 0;
		for (auto& p : ps)
		{
			sum += p.size();
		}
		sum /= ps.size();
		for (auto& p : ps)
		{
			if (totalLen < genomeSize)
			{
				generateContig(g, p, seq, seqOut, i++, totalLen);
			}
			std::copy(p.begin(), prev(p.end()), ostream_iterator<int, char>(outPath, " -> "));
			outPath << *p.rbegin() << endl;
		}
		//break;
		//set_union(un.begin(), un.end(), sp.begin(), sp.end(), inserter(un, un.begin()));
	}
	comps.clear();
	overlap.clear();
	outAssembly.close();
	outPath.close();
}

void cwd::readOverlapGraph(const string& graphName, vector<AGraph>& v_g)
{
	//using namespace boost;

	//// Construct an empty graph and prepare the dynamic_property_maps.
	//AGraph graph(0);
	//dynamic_properties dp(ignore_other_properties);

	//dp.property("label", get(&AVertex::r, graph));
	//dp.property("adj", get(&AEdge::adj, graph));
	//dp.property("weight", get(&AEdge::weight, graph));

	//// Sample graph as an std::istream;
	//ifstream dot(graphName);
	//FILE* fp = fopen(graphName.c_str(), "rb");
	//fread(&comps, 1, 1, fp);
	//std::istringstream
	//	gvgraph("digraph { graph [name=\"graphname\"]  a  c e [mass = 6.66] }");

	//while (read_graphviz(dot, graph, dp, "node_id"))
	//{
	//	v_g.push_back(graph);
	//}

}
