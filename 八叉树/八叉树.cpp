#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <cassert>
#include <set>
#include "RationalNumber.h"
#include "utility.h"

using namespace std;


struct OctTreeNode
{
	bool is_leaf;
	AABB box;
	vector<size_t> triangles;
	size_t num_of_tri_not_deleted = 0;
    vector<OctTreeNode*> child;
	OctTreeNode() :child(8, nullptr) {}
};



RationalNumber innerProduct(const PointOrVector2D& v1, const PointOrVector2D& v2)
{
	return v1.x * v2.x + v1.y * v2.y;
}

RationalNumber crossProduct(const PointOrVector2D& v1, const PointOrVector2D& v2) //大于等于0说明v1在v2右手，否则左手
{
	return v1.x * v2.y - v2.x * v1.y;
}

bool pointInConvexPolygon(const vector<PointOrVector2D>& h, const PointOrVector2D& p)
{
	for (size_t i = 0; i < h.size(); ++i)
	{
		PointOrVector2D pre = h[i];
		PointOrVector2D next = h[(i + 1) % h.size()];
		PointOrVector2D t = next - pre;
		PointOrVector2D n = p - pre;
		if (innerProduct(n, t) > r_zero)
		{
			return false;
		}
	}
	return true;
}

bool onOneLine(const vector<PointOrVector2D>& h)
{
	if (h.size() == 2)
		return true;
	for (size_t i = 0; i < h.size() - 2; ++i)
	{
		PointOrVector2D v1 = h[i + 1] - h[i];
		PointOrVector2D v2 = h[i + 2] - h[i];
		if (crossProduct(v1, v2) != r_zero)
		{
			return false;
		}
	}
	return true;
}

LineSegment getLineSegment(const vector<PointOrVector2D>& h)
{
	size_t min_index = 0;
	size_t max_index = 0;
	for (size_t i = 1; i < h.size(); ++i)
	{
		if (h[i].x < h[min_index].x)
		{
			min_index = i;
		}

		if (h[i].x > h[max_index].x)
		{
			max_index = i;
		}
	}

	if (h[min_index].x == h[max_index].x)
	{
		size_t min_index_y = 0;
		size_t max_index_y = 0;
		for (size_t i = 1; i < h.size(); ++i)
		{
			if (h[i].y < h[min_index_y].y)
			{
				min_index_y = i;
			}

			if (h[i].y > h[max_index_y].y)
			{
				max_index_y = i;
			}
		}
		return LineSegment(h[min_index_y].x, h[min_index_y].y, h[max_index_y].x, h[max_index_y].y);
	}
	else
	{
		return LineSegment(h[min_index].x, h[min_index].y, h[max_index].x, h[max_index].y);
	}
}

bool testPolygonInterForOne(const vector<PointOrVector2D>& left, const vector<PointOrVector2D>& right)
{
	for (size_t i = 0; i < left.size(); ++i)
	{
		PointOrVector2D pre = left[i];
		PointOrVector2D next = left[(i + 1) % left.size()];
		PointOrVector2D t = next - pre;
		PointOrVector2D n(-t.y, t.x);
		RationalNumber _min = innerProduct(left[0] - pre, n);
		RationalNumber _max = innerProduct(left[0] - pre, n);
		for (size_t j = 1; j < left.size(); ++j)
		{
			RationalNumber inner = innerProduct(left[j] - pre, n);
			_min = min(_min, inner);
			_max = max(_max, inner);
		}

		RationalNumber _minr = innerProduct(right[0] - pre, n);
		RationalNumber _maxr = innerProduct(right[0] - pre, n);
		for (size_t j = 1; j < right.size(); ++j)
		{
			RationalNumber inner = innerProduct(right[j] - pre, n);
			_minr = min(_minr, inner);
			_maxr = max(_maxr, inner);
		}

		if (_max < _minr || _maxr < _min)
		{
			return false;
		}
	}
	return true;
}

bool lineSegmentInterConvexPolygon(const LineSegment &l, const vector<PointOrVector2D>& h)
{
	if (pointInConvexPolygon(h, l.left) || pointInConvexPolygon(h, l.right))
	{
		return true;
	}

	for (size_t i = 0; i < h.size(); ++i)
	{
		PointOrVector2D pre = h[i];
		PointOrVector2D next = h[(i + 1) % h.size()];
		LineSegment edge(pre.x, pre.y, next.x, next.y);
		if (intersect(l, edge))
		{
			return true;
		}
	}
	return false;
}

bool ConvexPolygonInter(const vector<PointOrVector2D> &left, const vector<PointOrVector2D> &right)  //left,right均为凸多边形，顶点逆时针排列，该函数用分离轴定理判断它们是否相交
{
	if (onOneLine(left))
	{
		LineSegment l1 = getLineSegment(left);
		if (onOneLine(right))
		{
			LineSegment l2 = getLineSegment(right);
			return intersect(l1, l2);
		}
		return lineSegmentInterConvexPolygon(l1, right);
	}
	else
	{
		if (onOneLine(right))
		{
			LineSegment l2 = getLineSegment(right);
			return lineSegmentInterConvexPolygon(l2, left);
		}
		else
		{
			if (testPolygonInterForOne(left, right) == false)
			{
				return false;
			}
			return testPolygonInterForOne(right, left);
		}
	}
}

PointOrVector3D getProjPoint(PointOrVector3D &n, PointOrVector3D& d, PointOrVector3D& v, PointOrVector3D& p) //n要求交点的平面法向量,d投影直线方向向量,v平面上一点,p要投影的点,函数返回投影后的点
{

	PointOrVector3D s = subtract(v, p);
	RationalNumber l = innerProduct(n, s) / innerProduct(n, d);

	PointOrVector3D res;
	res.x = p.x + d.x * l;
	res.y = p.y + d.y * l;
	res.z = p.z + d.z * l;
	return res;
}

PointOrVector2D transformCoord(const PointOrVector3D& point, const PointOrVector3D& v, const PointOrVector3D& n, const PointOrVector3D _n)  //计算标架n,_n(彼此正交但不一定是单位向量)(在原直角坐标空间表示,原点为v)所在平面上的点point,在标价中的坐标
{
	PointOrVector3D _p_ = point - v;   //方程x * n + y * _n = point - v 注意n和_n为向量n= (n.x,n.y,n.z) _n(_n.x,_n.y,_n.z) _p_=(_p_.x, _p_.y,_p_.z) 方程展开后注意必有唯一解，找出系数矩阵线性无关的两行对应的方程，按克拉默法则求解
	RationalNumber nx;                 //这里可以point - v与n,_n做点积求得投影坐标，但这里采用了更复杂的做法,没有什么必要，仅供参考,修改为点积的做法当然很简单，可自己尝试
	RationalNumber ny;
	RationalNumber _nx;
	RationalNumber _ny;
	RationalNumber _px;
	RationalNumber _py;
	RationalNumber det;
	if ((det = n.x * _n.y - _n.x * n.y) != r_zero)
	{
		nx = n.x;
		ny = n.y;
		_nx = _n.x;
		_ny = _n.y;
		_px = _p_.x;
		_py = _p_.y;
	}
	else if ((det = n.x * _n.z - _n.x * n.z) != r_zero)
	{
		nx = n.x;
		ny = n.z;
		_nx = _n.x;
		_ny = _n.z;
		_px = _p_.x;
		_py = _p_.z;
	}
	else
	{
		nx = n.y;
		ny = n.z;
		_nx = _n.y;
		_ny = _n.z;
		_px = _p_.y;
		_py = _p_.z;
		det = nx * _ny - _nx * ny;
	}

	PointOrVector2D final_result_point;
	final_result_point.x = (_px * _ny - _nx * _py) / det;
	final_result_point.y = (nx * _py - _px * ny) / det;

	//final_result_point.x = innerProduct(_p_, n);
	//final_result_point.x = innerProduct(_p_, _n);  另一种计算方法,直接点积求得投影坐标,更简单,上面的方法仅供参考
	return  final_result_point;
}

bool equal(const RationalNumber& l, const RationalNumber& r)
{
	RationalNumber diff = l - r;
	diff.positive = true;
	if (diff < RationalNumber(1, 100000ll))
	{
		return true;
	}
	return false;
}

RationalNumber fastComputeSqrt(const RationalNumber& a)
{
	RationalNumber iter = a;
	while (equal(iter * iter, a) == false)
	{
		iter = (iter + a / iter) / RationalNumber(2, 1);
	}
	return iter;
}

bool triangleAndAABBIntersect(const AABB& box, const Triangle& t)  //判断三角形t和box是否相交
{
	PointOrVector3D r1 = t.v2 - t.v1;
	PointOrVector3D r2 = t.v3 - t.v1;
	PointOrVector3D n = crossProduct(r1, r2);

	PointOrVector3D s = crossProduct(n, r1);
	PointOrVector3D _n = crossProduct(s, n);

	RationalNumber len;
	RationalNumber _len;
	if ((_len = innerProduct(n, n)) != (len = innerProduct(_n, _n)))
	{
		RationalNumber len = fastComputeSqrt(len);
		RationalNumber _len = fastComputeSqrt(_len);
		_n.x = _n.x / len * _len;
		_n.y = _n.y / len * _len;
		_n.z = _n.z / len * _len;
	}

	vector<PointOrVector3D> AABBVertex { PointOrVector3D (box.xmax, box.ymax, box.zmax), PointOrVector3D(box.xmin, box.ymax, box.zmax), 
	PointOrVector3D(box.xmin, box.ymin, box.zmax) , PointOrVector3D(box.xmax, box.ymin, box.zmax) , PointOrVector3D(box.xmax, box.ymax, box.zmin) , 
		PointOrVector3D(box.xmin, box.ymax, box.zmin) , PointOrVector3D(box.xmin, box.ymin, box.zmin) , PointOrVector3D(box.xmax, box.ymin, box.zmin) };
	vector<PointOrVector3D> TriangleVertex { t.v1, t.v2, t.v3 };

	vector<PointOrVector3D> AABBProj(8);
	vector<PointOrVector3D> TriangleVertexProj(3);

	PointOrVector3D v1 = t.v1;
	for (size_t i = 0; i < 8; ++i)
	{
		AABBProj[i] = getProjPoint(s, s, v1, AABBVertex[i]);
	}

	for (size_t i = 0; i < 3; ++i)
	{
		TriangleVertexProj[i] = getProjPoint(s, s, v1, TriangleVertex[i]);
	}

	vector<PointOrVector2D> AABBProjOnPlanrCoord(8);
	vector<PointOrVector2D> TriangleVertexProjOnPlanrCoord(3); 

	for (size_t i = 0; i < 8; ++i)
	{
		AABBProjOnPlanrCoord[i] = transformCoord(AABBProj[i], v1, n, _n);
	}

	for (size_t i = 0; i < 3; ++i)
	{
		TriangleVertexProjOnPlanrCoord[i] = transformCoord(TriangleVertexProj[i], v1, n, _n);
	}

	vector<PointOrVector2D> ConvexPolygonAABB = getConvexHull(AABBProjOnPlanrCoord);
	vector<PointOrVector2D> ConvexPolygonTriangle = getConvexHull(TriangleVertexProjOnPlanrCoord);

	if (ConvexPolygonInter(ConvexPolygonAABB, ConvexPolygonTriangle) == false)  //AABB和三角形在与分离平面(三角形所在平面)垂直的投影平面上的投影凸包不相交,则三角形和AABB不相交
	{
		return false;
	}

	{    //三角形和AABB在xy,yz,xz中某一坐标平面上的投影不相交，则三角形和AABB不相交
		vector<PointOrVector2D> AABBProjxy{ PointOrVector2D(box.xmax, box.ymax), PointOrVector2D(box.xmin, box.ymax), PointOrVector2D(box.xmin, box.ymin), PointOrVector2D(box.xmax, box.ymin)};
		vector<PointOrVector2D> TriangleProjxy{ PointOrVector2D(t.v1.x, t.v1.y), PointOrVector2D(t.v2.x, t.v2.y), PointOrVector2D(t.v3.x, t.v3.y)};
		if (ConvexPolygonInter(AABBProjxy, TriangleProjxy) == false)
		{
			return false;
		}

		vector<PointOrVector2D> AABBProjxz{ PointOrVector2D(box.xmax, box.zmax), PointOrVector2D(box.xmin, box.zmax), PointOrVector2D(box.xmin, box.zmin), PointOrVector2D(box.xmax, box.zmin) };
		vector<PointOrVector2D> TriangleProjxz{ PointOrVector2D(t.v1.x, t.v1.z), PointOrVector2D(t.v2.x, t.v2.z), PointOrVector2D(t.v3.x, t.v3.z) };
		if (ConvexPolygonInter(AABBProjxz, TriangleProjxz) == false)
		{
			return false;
		}

		vector<PointOrVector2D> AABBProjyz{ PointOrVector2D(box.ymax, box.zmax), PointOrVector2D(box.ymin, box.zmax), PointOrVector2D(box.ymin, box.zmin), PointOrVector2D(box.ymax, box.zmin) };
		vector<PointOrVector2D> TriangleProjyz{ PointOrVector2D(t.v1.y, t.v1.z), PointOrVector2D(t.v2.y, t.v2.z), PointOrVector2D(t.v3.y, t.v3.z) };
		if (ConvexPolygonInter(AABBProjyz, TriangleProjyz) == false)
		{
			return false;
		}
	}
	return true;
}


int triangleAndAABBIntersectPrecise(const AABB& box, const Triangle& t)
{
	if (triangleAndAABBIntersect(box, t) == false)   //判断三角形是否完全位于AABB盒外,如果完全位于AABB盒外返回1,否则继续
		return 1;
	RationalNumber _xmin = min({ t.v1.x, t.v2.x, t.v3.x });
	RationalNumber _xmax = max({ t.v1.x, t.v2.x, t.v3.x });
	RationalNumber _ymin = min({ t.v1.y, t.v2.y, t.v3.y });
	RationalNumber _ymax = max({ t.v1.y, t.v2.y, t.v3.y });
	RationalNumber _zmin = min({ t.v1.z, t.v2.z, t.v3.z });
	RationalNumber _zmax = max({ t.v1.z, t.v2.z, t.v3.z });

	if (box.xmin <= _xmin && _xmax <= box.xmax && box.ymin <= _ymin && _ymax <= box.ymax && box.zmin <= _zmin && _zmax <= box.zmax) //三角形位于AABB盒内(三角形位于AABB内当且仅当其包围盒位于AABB内)
	{
		return 2;
	}
	return 3;   //t一部分在box外部,一部分完全位于box内部或一部分位于box外部,另一部分在box表面
}

bool triangleHasPartStrictInAABB(const AABB& box, const Triangle& t)  //注意,这里假定box和t相交,判断t是否交于box内部,只要三角形与AABB相交,则三角形的包围盒在xyz三条坐标轴上的投影区间和AABB对应投影区间均不只交于一点,则三角形
{                                                          //交于AABB内部,否则就不是
	RationalNumber _xmin = min({ t.v1.x, t.v2.x, t.v3.x });   
	RationalNumber _xmax = max({ t.v1.x, t.v2.x, t.v3.x });
	RationalNumber _ymin = min({ t.v1.y, t.v2.y, t.v3.y });
	RationalNumber _ymax = max({ t.v1.y, t.v2.y, t.v3.y });
	RationalNumber _zmin = min({ t.v1.z, t.v2.z, t.v3.z });
	RationalNumber _zmax = max({ t.v1.z, t.v2.z, t.v3.z });

	bool inter;
	pair<Value, Value> r = getIntersection(pair<Value, Value>(Value(_xmin), Value(_xmax)), pair<Value, Value>(Value(box.xmin), Value(box.xmax)), inter);

	if (r.first == r.second)
	{
		return false;
	}

	r = getIntersection(pair<Value, Value>(Value(_ymin), Value(_ymax)), pair<Value, Value>(Value(box.ymin), Value(box.ymax)), inter);


	if (r.first == r.second)
	{
		return false;
	}

	r = getIntersection(pair<Value, Value>(Value(_zmin), Value(_zmax)), pair<Value, Value>(Value(box.zmin), Value(box.zmax)), inter);

	if (r.first == r.second)
	{
		return false;
	}
	return true;
}

enum InSurface { XY, XZ, YZ };
bool triangleInterAABBWithOneSurface(const AABB& box, const Triangle& t, InSurface s)  //这里假定box,t相交，并且t和box某个面共面,判断t是否只和box某一面交于二维面片
{
	switch (s)     //这里逻辑和triangleHasPartStrictInAABB函数类似，只不过需要把t和t共面的box面投影到对应坐标面再判断投影是否交于内部(转换为二维问题),s是需要投影到的的坐标面
	{
	case InSurface::XY:
	{
		RationalNumber _xmin = min({ t.v1.x, t.v2.x, t.v3.x });
		RationalNumber _xmax = max({ t.v1.x, t.v2.x, t.v3.x });
		RationalNumber _ymin = min({ t.v1.y, t.v2.y, t.v3.y });
		RationalNumber _ymax = max({ t.v1.y, t.v2.y, t.v3.y });

		bool inter;
		pair<Value, Value> r = getIntersection(pair<Value, Value>(Value(_xmin), Value(_xmax)), pair<Value, Value>(Value(box.xmin), Value(box.xmax)), inter);

		if (r.first == r.second)
		{
			return false;
		}

		r = getIntersection(pair<Value, Value>(Value(_ymin), Value(_ymax)), pair<Value, Value>(Value(box.ymin), Value(box.ymax)), inter);


		if (r.first == r.second)
		{
			return false;
		}
		break;
	}
	case InSurface::XZ:
	{
		RationalNumber _xmin = min({ t.v1.x, t.v2.x, t.v3.x });
		RationalNumber _xmax = max({ t.v1.x, t.v2.x, t.v3.x });
		RationalNumber _zmin = min({ t.v1.z, t.v2.z, t.v3.z });
		RationalNumber _zmax = max({ t.v1.z, t.v2.z, t.v3.z });

		bool inter;
		pair<Value, Value> r = getIntersection(pair<Value, Value>(Value(_xmin), Value(_xmax)), pair<Value, Value>(Value(box.xmin), Value(box.xmax)), inter);

		if (r.first == r.second)
		{
			return false;
		}

		r = getIntersection(pair<Value, Value>(Value(_zmin), Value(_zmax)), pair<Value, Value>(Value(box.zmin), Value(box.zmax)), inter);


		if (r.first == r.second)
		{
			return false;
		}
		break;
	}
	case InSurface::YZ:
	{
		RationalNumber _ymin = min({ t.v1.y, t.v2.y, t.v3.y });
		RationalNumber _ymax = max({ t.v1.y, t.v2.y, t.v3.y });
		RationalNumber _zmin = min({ t.v1.z, t.v2.z, t.v3.z });
		RationalNumber _zmax = max({ t.v1.z, t.v2.z, t.v3.z });

		bool inter;
		pair<Value, Value> r = getIntersection(pair<Value, Value>(Value(_ymin), Value(_ymax)), pair<Value, Value>(Value(box.ymin), Value(box.ymax)), inter);

		if (r.first == r.second)
		{
			return false;
		}

		r = getIntersection(pair<Value, Value>(Value(_zmin), Value(_zmax)), pair<Value, Value>(Value(box.zmin), Value(box.zmax)), inter);


		if (r.first == r.second)
		{
			return false;
		}
		break;
	}
	}
	return true;
}

class OctTree
{
public:
	OctTree(const vector<Triangle>& triangles, size_t r, size_t c);
	OctTree(size_t r, size_t c) :recursive_max_depth(r), capacity_min(c), triangles(), is_removable(){}
	pair<std::vector<pair<size_t, InterType>>, Value> solveLineAndTriangleInter(PointOrVector3D& p0, PointOrVector3D& d, Value tmin, Value tmax);
	bool removeTriangle(size_t index);
	bool insert(const Triangle& t);
	~OctTree()
	{
		deleteOctTree(root);
	}
	const vector<Triangle> &getTriangles() const
	{
		return triangles;
	}
private:
    void buildOctTree(const vector<size_t>& input, OctTreeNode* root, size_t depth, bool insert_mode_first_enter);
	void insertTriangleToOctTree(OctTreeNode* root, const Triangle& t, bool insert_be_deleted, size_t index);
	bool remove(OctTreeNode *cur, size_t index);
	void insertIntoSubTree(OctTreeNode* cur, size_t i, const AABB& b, size_t index, const Triangle& t, bool insert_be_deleted);
	void removeInSubTree(OctTreeNode* cur, size_t i, size_t index, bool& all_null);
	void deleteOctTree(OctTreeNode* root)
	{
		if (root == nullptr)
			return;
		for (size_t i = 0; i < 8; ++i)
		{
			deleteOctTree(root->child[i]);
		}
		delete root;
	}
	pair<std::vector<pair<size_t, InterType>>, Value> searchInter(OctTreeNode* cur, PointOrVector3D& p0, PointOrVector3D& d, Value tmin, Value tmax);
	OctTreeNode* root = nullptr;
    size_t recursive_max_depth;
	size_t capacity_min;
	vector<Triangle> triangles;
	vector<int> is_removable;
};

void OctTree::removeInSubTree(OctTreeNode* cur, size_t i, size_t index, bool &all_null)
{
	bool _r = remove(cur->child[i], index);

	if (_r == false)
	{
		all_null = false;
	}
	else
	{
		cur->child[i] = nullptr;
	}
}

bool OctTree::remove(OctTreeNode* cur, size_t index) //调用时cur总不为空,且被删除三角形总和cur相交,true表示cur节点已删除,false未删除
{
	if (cur->is_leaf)
	{
		--(cur->num_of_tri_not_deleted);
		if (cur->num_of_tri_not_deleted == 0)
		{
			delete cur;
			return true;
		}
		return false;
	}

	bool all_null = true;
	RationalNumber xmid = (cur->box.xmin + cur->box.xmax) / RationalNumber(2, 1);
	RationalNumber ymid = (cur->box.ymin + cur->box.ymax) / RationalNumber(2, 1);
	RationalNumber zmid = (cur->box.zmin + cur->box.zmax) / RationalNumber(2, 1);
	vector<int> tri_exclude(8, 0);

	for (size_t i = 0; i < 8; ++i)
	{
		if (cur->child[i] == nullptr)
			continue;

		if (tri_exclude[i])
		{
			all_null = false;
			continue;
		}
			

		AABB b = cur->child[i]->box;

		switch (triangleAndAABBIntersectPrecise(b, triangles[index]))
		{
		case 1:
		{
			all_null = false;
			break;
		}
		case 2:
		{
			for (size_t m = i + 1; m < 8; ++m)
			{
				tri_exclude[m] = 1;
			}
			removeInSubTree(cur, i, index, all_null);
			break;
		}
		case 3:
		{
			if (triangleHasPartStrictInAABB(b, triangles[index]))
			{
				removeInSubTree(cur, i, index, all_null);
			}
			else
			{
				switch (i)
				{
				case 0:
				{
					if (triangles[index].v1.x == xmid && triangles[index].v2.x == xmid && triangles[index].v3.x == xmid)
					{
						tri_exclude[1] = 1;
						removeInSubTree(cur, i, index, all_null);
						continue;
					}
					else if (triangles[index].v1.y == ymid && triangles[index].v2.y == ymid && triangles[index].v3.y == ymid)
					{
						tri_exclude[3] = 1;
						removeInSubTree(cur, i, index, all_null);
						continue;
					}
					else if (triangles[index].v1.z == zmid && triangles[index].v2.z == zmid && triangles[index].v3.z == zmid)
					{
						tri_exclude[4] = 1;
						removeInSubTree(cur, i, index, all_null);
						continue;
					}
					else if (triangles[index].v1.x == cur->box.xmax && triangles[index].v2.x == cur->box.xmax && triangles[index].v3.x == cur->box.xmax
						|| triangles[index].v1.y == cur->box.ymax && triangles[index].v2.y == cur->box.ymax && triangles[index].v3.y == cur->box.ymax
						|| triangles[index].v1.z == cur->box.zmax && triangles[index].v2.z == cur->box.zmax && triangles[index].v3.z == cur->box.zmax)
					{
						if (triangles[index].v1.x == cur->box.xmax && triangles[index].v2.x == cur->box.xmax && triangles[index].v3.x == cur->box.xmax)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::YZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else if (triangles[index].v1.y == cur->box.ymax && triangles[index].v2.y == cur->box.ymax && triangles[index].v3.y == cur->box.ymax)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XY))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
					}
					break;
				}
				case 1:
				{
					if (triangles[index].v1.y == ymid && triangles[index].v2.y == ymid && triangles[index].v3.y == ymid)
					{
						tri_exclude[2] = 1;
						removeInSubTree(cur, i, index, all_null);
					}
					else if (triangles[index].v1.z == zmid && triangles[index].v2.z == zmid && triangles[index].v3.z == zmid)
					{
						tri_exclude[5] = 1;
						removeInSubTree(cur, i, index, all_null);
					}
					else if (triangles[index].v1.x == cur->box.xmin && triangles[index].v2.x == cur->box.xmin && triangles[index].v3.x == cur->box.xmin
						|| triangles[index].v1.y == cur->box.ymax && triangles[index].v2.y == cur->box.ymax && triangles[index].v3.y == cur->box.ymax
						|| triangles[index].v1.z == cur->box.zmax && triangles[index].v2.z == cur->box.zmax && triangles[index].v3.z == cur->box.zmax)
					{
						if (triangles[index].v1.x == cur->box.xmin && triangles[index].v2.x == cur->box.xmin && triangles[index].v3.x == cur->box.xmin)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::YZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else if (triangles[index].v1.y == cur->box.ymax && triangles[index].v2.y == cur->box.ymax && triangles[index].v3.y == cur->box.ymax)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XY))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
					}
					break;
				}
				case 2:
				{
					if (triangles[index].v1.x == xmid && triangles[index].v2.x == xmid && triangles[index].v3.x == xmid)
					{
						tri_exclude[3] = 1;
						removeInSubTree(cur, i, index, all_null);
					}
					else if (triangles[index].v1.z == zmid && triangles[index].v2.z == zmid && triangles[index].v3.z == zmid)
					{
						tri_exclude[6] = 1;
						removeInSubTree(cur, i, index, all_null);
					}
					else if (triangles[index].v1.x == cur->box.xmin && triangles[index].v2.x == cur->box.xmin && triangles[index].v3.x == cur->box.xmin
						|| triangles[index].v1.y == cur->box.ymin && triangles[index].v2.y == cur->box.ymin && triangles[index].v3.y == cur->box.ymin
						|| triangles[index].v1.z == cur->box.zmax && triangles[index].v2.z == cur->box.zmax && triangles[index].v3.z == cur->box.zmax)
					{
						if (triangles[index].v1.x == cur->box.xmin && triangles[index].v2.x == cur->box.xmin && triangles[index].v3.x == cur->box.xmin)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::YZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else if (triangles[index].v1.y == cur->box.ymin && triangles[index].v2.y == cur->box.ymin && triangles[index].v3.y == cur->box.ymin)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XY))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
					}
					break;
				}
				case 3:
				{
					if (triangles[index].v1.z == zmid && triangles[index].v2.z == zmid && triangles[index].v3.z == zmid)
					{
						tri_exclude[7] = 1;
						removeInSubTree(cur, i, index, all_null);
					}
					else if (triangles[index].v1.x == cur->box.xmax && triangles[index].v2.x == cur->box.xmax && triangles[index].v3.x == cur->box.xmax
						|| triangles[index].v1.y == cur->box.ymin && triangles[index].v2.y == cur->box.ymin && triangles[index].v3.y == cur->box.ymin
						|| triangles[index].v1.z == cur->box.zmax && triangles[index].v2.z == cur->box.zmax && triangles[index].v3.z == cur->box.zmax)
					{
						if (triangles[index].v1.x == cur->box.xmax && triangles[index].v2.x == cur->box.xmax && triangles[index].v3.x == cur->box.xmax)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::YZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else if (triangles[index].v1.y == cur->box.ymin && triangles[index].v2.y == cur->box.ymin && triangles[index].v3.y == cur->box.ymin)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XY))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
					}
					break;
				}
				case 4:
				{
					if (triangles[index].v1.x == xmid && triangles[index].v2.x == xmid && triangles[index].v3.x == xmid)
					{
						tri_exclude[5] = 1;
						removeInSubTree(cur, i, index, all_null);
					}
					else if (triangles[index].v1.y == ymid && triangles[index].v2.y == ymid && triangles[index].v3.y == ymid)
					{
						tri_exclude[7] = 1;
						removeInSubTree(cur, i, index, all_null);
					}
					else if (triangles[index].v1.x == cur->box.xmax && triangles[index].v2.x == cur->box.xmax && triangles[index].v3.x == cur->box.xmax
						|| triangles[index].v1.y == cur->box.ymax && triangles[index].v2.y == cur->box.ymax && triangles[index].v3.y == cur->box.ymax
						|| triangles[index].v1.z == cur->box.zmin && triangles[index].v2.z == cur->box.zmin && triangles[index].v3.z == cur->box.zmin)
					{
						if (triangles[index].v1.x == cur->box.xmax && triangles[index].v2.x == cur->box.xmax && triangles[index].v3.x == cur->box.xmax)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::YZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else if (triangles[index].v1.y == cur->box.ymax && triangles[index].v2.y == cur->box.ymax && triangles[index].v3.y == cur->box.ymax)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XY))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
					}
					break;
				}
				case 5:
				{
					if (triangles[index].v1.y == ymid && triangles[index].v2.y == ymid && triangles[index].v3.y == ymid)
					{
						tri_exclude[6] = 1;
						removeInSubTree(cur, i, index, all_null);
					}
					else if (triangles[index].v1.x == cur->box.xmin && triangles[index].v2.x == cur->box.xmin && triangles[index].v3.x == cur->box.xmin
						|| triangles[index].v1.y == cur->box.ymax && triangles[index].v2.y == cur->box.ymax && triangles[index].v3.y == cur->box.ymax
						|| triangles[index].v1.z == cur->box.zmin && triangles[index].v2.z == cur->box.zmin && triangles[index].v3.z == cur->box.zmin)
					{
						if (triangles[index].v1.x == cur->box.xmin && triangles[index].v2.x == cur->box.xmin && triangles[index].v3.x == cur->box.xmin)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::YZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else if (triangles[index].v1.y == cur->box.ymax && triangles[index].v2.y == cur->box.ymax && triangles[index].v3.y == cur->box.ymax)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XY))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
					}
					break;
				}
				case 6:
				{
					if (triangles[index].v1.x == xmid && triangles[index].v2.x == xmid && triangles[index].v3.x == xmid)
					{
						tri_exclude[7] = 1;
						removeInSubTree(cur, i, index, all_null);
					}
					else if (triangles[index].v1.x == cur->box.xmin && triangles[index].v2.x == cur->box.xmin && triangles[index].v3.x == cur->box.xmin
						|| triangles[index].v1.y == cur->box.ymin && triangles[index].v2.y == cur->box.ymin && triangles[index].v3.y == cur->box.ymin
						|| triangles[index].v1.z == cur->box.zmin && triangles[index].v2.z == cur->box.zmin && triangles[index].v3.z == cur->box.zmin)
					{
						if (triangles[index].v1.x == cur->box.xmin && triangles[index].v2.x == cur->box.xmin && triangles[index].v3.x == cur->box.xmin)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::YZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else if (triangles[index].v1.y == cur->box.ymin && triangles[index].v2.y == cur->box.ymin && triangles[index].v3.y == cur->box.ymin)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XY))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
					}
					break;
				}
				case 7:
				{
					if (triangles[index].v1.x == cur->box.xmax && triangles[index].v2.x == cur->box.xmax && triangles[index].v3.x == cur->box.xmax
						|| triangles[index].v1.y == cur->box.ymin && triangles[index].v2.y == cur->box.ymin && triangles[index].v3.y == cur->box.ymin
						|| triangles[index].v1.z == cur->box.zmin && triangles[index].v2.z == cur->box.zmin && triangles[index].v3.z == cur->box.zmin)
					{
						if (triangles[index].v1.x == cur->box.xmax && triangles[index].v2.x == cur->box.xmax && triangles[index].v3.x == cur->box.xmax)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::YZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else if (triangles[index].v1.y == cur->box.ymin && triangles[index].v2.y == cur->box.ymin && triangles[index].v3.y == cur->box.ymin)
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XZ))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
						else
						{
							if (triangleInterAABBWithOneSurface(b, triangles[index], InSurface::XY))
							{
								removeInSubTree(cur, i, index, all_null);
							}
						}
					}
					break;
				}
				}
			}

			break;
		}
		}
	}

	if (all_null)
	{
		delete cur;
		return true;
	}
	return false;
}

bool OctTree::removeTriangle(size_t index)
{
	if (index >= triangles.size())
		return false;
	if (is_removable[index] == 1)
		return false;
	is_removable[index] = 1;

	bool _r = remove(root, index);
	if (_r)
		root = nullptr;
	return true;
}

//Value是最近相交点的p+t*d的t值,vector<pair(存放和线段相交的三角形所在triangle数组的索引,和三角形的相交类型)>
pair<std::vector<pair<size_t, InterType>>, Value> OctTree::searchInter(OctTreeNode* cur, PointOrVector3D& p0, PointOrVector3D& d, Value tmin, Value tmax)  //要求tmin<=tmax
{
	pair<Value, Value> res = judgeAABBInterline(cur->box.xmin, cur->box.ymin, cur->box.zmin, cur->box.xmax, cur->box.ymax, cur->box.zmax,
		p0.x, p0.y, p0.z, d.x, d.y, d.z, tmin.value, tmax.value, tmin.type, tmax.type);
	if (res.first == v_type::positive_inf)
	{
	    return { std::vector<pair<size_t, InterType>>(), Value(v_type::positive_inf) };
	}

	if (cur->is_leaf)   //这里可以把叶节点的三角形放入BVH树并搜索BVH树,这里为简单起见没有这样做
	{
		std::vector<pair<size_t, InterType>> _r;
		_r.reserve(cur->triangles.size());
		Value t(v_type::positive_inf);
		for (size_t i = 0; i < cur->triangles.size(); ++i)
		{
			if (is_removable[cur->triangles[i]] == 1)
			    continue;
			pair<RationalNumber, InterType> _res = lineInterTriangle(p0, d, tmin, tmax, triangles[cur->triangles[i]].v1, triangles[cur->triangles[i]].v2, triangles[cur->triangles[i]].v3);
			if (_res.second != InterType::NOInter)
			{
				if (Value(_res.first) < t)
				{
					t = Value(_res.first);
					_r.clear();
					_r.push_back({ cur->triangles[i], _res.second });
				}
				else if (Value(_res.first) == t)
				{
					_r.push_back({ cur->triangles[i], _res.second });
				}
			}
		}
		_r.shrink_to_fit();
		return { _r , t };
	}

	pair<std::vector<pair<size_t, InterType>>, Value> result;
	result.second = Value(v_type::positive_inf);
	for (size_t i = 0; i < 8; ++i)
	{
		if (cur->child[i] == nullptr)
			continue;
		pair<std::vector<pair<size_t, InterType>>, Value> _r = searchInter(cur->child[i], p0, d, res.first, res.second);

		if (_r.second.type == v_type::positive_inf)
			continue;

		if (_r.second < result.second)
		{
			result = _r;
		}
		else if (_r.second == result.second)
		{
			vector<int> work(triangles.size(), 0);
			for (size_t j = 0; j < result.first.size(); ++j)
			{
				work[result.first[j].first] = 1;
			}

			for (size_t j = 0; j < _r.first.size(); ++j)
			{
				if (work[_r.first[j].first] == 0)
				{
					result.first.push_back(_r.first[j]);
				}
			}
		}
	}
	return result;
}

pair<std::vector<pair<size_t, InterType>>, Value> OctTree::solveLineAndTriangleInter(PointOrVector3D& p0, PointOrVector3D& d, Value tmin, Value tmax)
{
	if (root == nullptr)
	{
		return { std::vector<pair<size_t, InterType>>(), Value(v_type::positive_inf) };
	}
	pair<std::vector<pair<size_t, InterType>>, Value> res = searchInter(root, p0, d, tmin, tmax);
	return res;
}

OctTree::OctTree(const vector<Triangle>& triangles, size_t r, size_t c) :triangles(triangles), recursive_max_depth(r), capacity_min(c), is_removable(triangles.size())  //triangles不能为空不能重复,且r必须大于0, capacity_min必须大于0
{
	RationalNumber xmin = min({ triangles[0].v1.x, triangles[0].v2.x, triangles[0].v3.x });
	RationalNumber ymin = min({ triangles[0].v1.y, triangles[0].v2.y, triangles[0].v3.y });
	RationalNumber zmin = min({ triangles[0].v1.z, triangles[0].v2.z, triangles[0].v3.z });
	RationalNumber xmax = max({ triangles[0].v1.x, triangles[0].v2.x, triangles[0].v3.x });
	RationalNumber ymax = max({ triangles[0].v1.y, triangles[0].v2.y, triangles[0].v3.y });
	RationalNumber zmax = max({ triangles[0].v1.z, triangles[0].v2.z, triangles[0].v3.z });
	for (size_t i = 1; i < triangles.size(); ++i)
	{
		xmin = min({ triangles[i].v1.x, triangles[i].v2.x, triangles[i].v3.x, xmin });
		ymin = min({ triangles[i].v1.y, triangles[i].v2.y, triangles[i].v3.y, ymin });
		zmin = min({ triangles[i].v1.z, triangles[i].v2.z, triangles[i].v3.z, zmin });
		xmax = max({ triangles[i].v1.x, triangles[i].v2.x, triangles[i].v3.x, xmax });
		ymax = max({ triangles[i].v1.y, triangles[i].v2.y, triangles[i].v3.y, ymax });
		zmax = max({ triangles[i].v1.z, triangles[i].v2.z, triangles[i].v3.z, zmax });
	}

	AABB b_box;
	b_box.xmin = xmin;
	b_box.ymin = ymin;
	b_box.zmin = zmin;
	b_box.xmax = xmax;
	b_box.ymax = ymax;
	b_box.zmax = zmax;

	root = new OctTreeNode();
	root->box = b_box;

	vector<size_t> input(triangles.size());
	for (size_t i = 0; i < input.size(); ++i)
	{
		input[i] = i;
	}
	buildOctTree(input, root, 1, false);
}

void computeAABB(AABB &_box, RationalNumber xmid, RationalNumber ymid, RationalNumber zmid, const AABB &p, size_t i)
{
	switch (i)   //sub[0]_max p            sub[0]_min mid             root_min p
	{
	    case 0:
		{
			_box.xmin = xmid;
			_box.ymin = ymid;
			_box.zmin = zmid;

			_box.xmax = p.xmax;
			_box.ymax = p.ymax;
			_box.zmax = p.zmax;
			break;
		}
		case 1:
		{
			_box.xmin = p.xmin;
			_box.ymin = ymid;
			_box.zmin = zmid;

			_box.xmax = xmid;
			_box.ymax = p.ymax;
			_box.zmax = p.zmax;
			break;
		}
		case 2:
		{
			_box.xmin = p.xmin;
			_box.ymin = p.ymin;
			_box.zmin = zmid;

			_box.xmax = xmid;
			_box.ymax = ymid;
			_box.zmax = p.zmax;
			break;
		}
		case 3:
		{
			_box.xmin = xmid;
			_box.ymin = p.ymin;
			_box.zmin = zmid;

			_box.xmax = p.xmax;
			_box.ymax = ymid;
			_box.zmax = p.zmax;
			break;
		}
		case 4:
		{
			_box.xmin = xmid;
			_box.ymin = ymid;
			_box.zmin = p.zmin;

			_box.xmax = p.xmax;
			_box.ymax = p.ymax;
			_box.zmax = zmid;
			break;
		}
		case 5:
		{
			_box.xmin = p.xmin;
			_box.ymin = ymid;
			_box.zmin = p.zmin;

			_box.xmax = xmid;
			_box.ymax = p.ymax;
			_box.zmax = zmid;
			break;
		}
		case 6:
		{
			_box.xmin = p.xmin;
			_box.ymin = p.ymin;
			_box.zmin = p.zmin;

			_box.xmax = xmid;
			_box.ymax = ymid;
			_box.zmax = zmid;
			break;
		}
		case 7:
		{
			_box.xmin = xmid;
			_box.ymin = p.ymin;
			_box.zmin = p.zmin;

			_box.xmax = p.xmax;
			_box.ymax = ymid;
			_box.zmax = zmid;
			break;
		}
	}

}

void OctTree::buildOctTree(const vector<size_t>& input, OctTreeNode* root, size_t depth, bool insert_mode_first_enter)
{
	if ((insert_mode_first_enter ? root->num_of_tri_not_deleted : input.size()) <= capacity_min || depth == recursive_max_depth)
	{
		root->is_leaf = true;
		root->triangles = input;
		if (root->num_of_tri_not_deleted == 0)
		{
			root->num_of_tri_not_deleted = input.size();
		}
		return;
	}

	root->is_leaf = false;
	RationalNumber xmid = (root->box.xmin + root->box.xmax) / RationalNumber(2, 1);
	RationalNumber ymid = (root->box.ymin + root->box.ymax) / RationalNumber(2, 1);
	RationalNumber zmid = (root->box.zmin + root->box.zmax) / RationalNumber(2, 1);
	vector<AABB> sub_box(8);  //依次为+x+y+z -x+y+z -x-y+z +x-y+z +x+y-z -x+y-z -x-y-z +x-y-z
	computeAABB(sub_box[0], xmid, ymid, zmid, root->box, 0);
	computeAABB(sub_box[1], xmid, ymid, zmid, root->box, 1);
	computeAABB(sub_box[2], xmid, ymid, zmid, root->box, 2);
	computeAABB(sub_box[3], xmid, ymid, zmid, root->box, 3);
	computeAABB(sub_box[4], xmid, ymid, zmid, root->box, 4);
	computeAABB(sub_box[5], xmid, ymid, zmid, root->box, 5);
	computeAABB(sub_box[6], xmid, ymid, zmid, root->box, 6);
	computeAABB(sub_box[7], xmid, ymid, zmid, root->box, 7);

	vector<vector<size_t>> input_t(8, vector<size_t>());
	vector<vector<int>> tri_exclude(input.size(), vector<int>(8, 0));
	for (size_t i = 0; i < 8; ++i)
	{
		input_t[i].reserve(root->num_of_tri_not_deleted == 0 ? input.size() : root->num_of_tri_not_deleted);
		for (size_t k = 0; k < input.size(); ++k)
		{
			if (is_removable[input[k]] == 1)
				continue;
			if (tri_exclude[k][i])  //索引input[k]对应三角形被小于i的子包围盒排除,不应该被分配到第i个子包围盒
				continue;
			switch (triangleAndAABBIntersectPrecise(sub_box[i], triangles[input[k]]))
			{
			case 1:
			{
				break;  //input[k]不分配给第i个子包围盒,i也不排除编号更大的子包围盒
			}
			case 2:
			{
				for (size_t m = i + 1; m < 8; ++m)//input[k]分配给第i个子包围盒,i排除所有编号更大的子包围盒
				{
					tri_exclude[k][m] = 1;
				}
				input_t[i].push_back(input[k]);
				break;
			}
			case 3:   // 
			{   

				if (triangleHasPartStrictInAABB(sub_box[i], triangles[input[k]]))  //input[k]一部分在box外部, 一部分完全位于box内部,此时将三角形分配给box
				{
					input_t[i].push_back(input[k]);
				}
				else
				{
					switch (i)  //最后一个else if之前的判断处理三角形input[k]在box的三个内部面所在平面上且和内部面交于平面，线段或点的情形
					{  //最后一个else if处理三角形input[k]在box的三个外部面所在平面上且和外部面交于平面，线段或点的情形(交于面则分配input[k]给box,否则什么都不做),当然剩下的情形是input[k]不在box六个面所在平面上且和box仅交于点线结构
					case 0:   //此时input[k]必有一顶点完全位于box外,对于这种情形我们什么都不做
					{
						if (triangles[input[k]].v1.x == xmid && triangles[input[k]].v2.x == xmid && triangles[input[k]].v3.x == xmid)
						{
							tri_exclude[k][1] = 1;
							input_t[i].push_back(input[k]);
							continue;
						}
						else if (triangles[input[k]].v1.y == ymid && triangles[input[k]].v2.y == ymid && triangles[input[k]].v3.y == ymid)
						{
							tri_exclude[k][3] = 1;
							input_t[i].push_back(input[k]);
							continue;
						}
						else if (triangles[input[k]].v1.z == zmid && triangles[input[k]].v2.z == zmid && triangles[input[k]].v3.z == zmid)
						{
							tri_exclude[k][4] = 1;
							input_t[i].push_back(input[k]);
							continue;
						}
						else if (triangles[input[k]].v1.x == root->box.xmax && triangles[input[k]].v2.x == root->box.xmax && triangles[input[k]].v3.x == root->box.xmax
							|| triangles[input[k]].v1.y == root->box.ymax && triangles[input[k]].v2.y == root->box.ymax && triangles[input[k]].v3.y == root->box.ymax
							|| triangles[input[k]].v1.z == root->box.zmax && triangles[input[k]].v2.z == root->box.zmax && triangles[input[k]].v3.z == root->box.zmax)
						{
							if (triangles[input[k]].v1.x == root->box.xmax && triangles[input[k]].v2.x == root->box.xmax && triangles[input[k]].v3.x == root->box.xmax)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::YZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else if (triangles[input[k]].v1.y == root->box.ymax && triangles[input[k]].v2.y == root->box.ymax && triangles[input[k]].v3.y == root->box.ymax)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XY))
								{
									input_t[i].push_back(input[k]);
								}
							}
						}
						break;
					}
					case 1:  //这里及之后没有写出的if else if 判断(不考虑的box的内部面)是不会通过的,因为input[k]已经被小于i的子包围盒排除
					{
						if (triangles[input[k]].v1.y == ymid && triangles[input[k]].v2.y == ymid && triangles[input[k]].v3.y == ymid)
						{
							tri_exclude[k][2] = 1;
							input_t[i].push_back(input[k]);
						}
						else if (triangles[input[k]].v1.z == zmid && triangles[input[k]].v2.z == zmid && triangles[input[k]].v3.z == zmid)
						{
							tri_exclude[k][5] = 1;
							input_t[i].push_back(input[k]);
						}
						else if (triangles[input[k]].v1.x == root->box.xmin && triangles[input[k]].v2.x == root->box.xmin && triangles[input[k]].v3.x == root->box.xmin
							|| triangles[input[k]].v1.y == root->box.ymax && triangles[input[k]].v2.y == root->box.ymax && triangles[input[k]].v3.y == root->box.ymax
							|| triangles[input[k]].v1.z == root->box.zmax && triangles[input[k]].v2.z == root->box.zmax && triangles[input[k]].v3.z == root->box.zmax)
						{
							if (triangles[input[k]].v1.x == root->box.xmin && triangles[input[k]].v2.x == root->box.xmin && triangles[input[k]].v3.x == root->box.xmin)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::YZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else if (triangles[input[k]].v1.y == root->box.ymax && triangles[input[k]].v2.y == root->box.ymax && triangles[input[k]].v3.y == root->box.ymax)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XY))
								{
									input_t[i].push_back(input[k]);
								}
							}
						}
						break;
					}
					case 2:
					{
						if (triangles[input[k]].v1.x == xmid && triangles[input[k]].v2.x == xmid && triangles[input[k]].v3.x == xmid)
						{
							tri_exclude[k][3] = 1;
							input_t[i].push_back(input[k]);
						}
						else if (triangles[input[k]].v1.z == zmid && triangles[input[k]].v2.z == zmid && triangles[input[k]].v3.z == zmid)
						{
							tri_exclude[k][6] = 1;
							input_t[i].push_back(input[k]);
						}
						else if (triangles[input[k]].v1.x == root->box.xmin && triangles[input[k]].v2.x == root->box.xmin && triangles[input[k]].v3.x == root->box.xmin
							|| triangles[input[k]].v1.y == root->box.ymin && triangles[input[k]].v2.y == root->box.ymin && triangles[input[k]].v3.y == root->box.ymin
							|| triangles[input[k]].v1.z == root->box.zmax && triangles[input[k]].v2.z == root->box.zmax && triangles[input[k]].v3.z == root->box.zmax)
						{
							if (triangles[input[k]].v1.x == root->box.xmin && triangles[input[k]].v2.x == root->box.xmin && triangles[input[k]].v3.x == root->box.xmin)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::YZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else if (triangles[input[k]].v1.y == root->box.ymin && triangles[input[k]].v2.y == root->box.ymin && triangles[input[k]].v3.y == root->box.ymin)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XY))
								{
									input_t[i].push_back(input[k]);
								}
							}
						}
						break;
					}
					case 3:
					{
						if (triangles[input[k]].v1.z == zmid && triangles[input[k]].v2.z == zmid && triangles[input[k]].v3.z == zmid)
						{
							tri_exclude[k][7] = 1;
							input_t[i].push_back(input[k]);
						}
						else if (triangles[input[k]].v1.x == root->box.xmax && triangles[input[k]].v2.x == root->box.xmax && triangles[input[k]].v3.x == root->box.xmax
							|| triangles[input[k]].v1.y == root->box.ymin && triangles[input[k]].v2.y == root->box.ymin && triangles[input[k]].v3.y == root->box.ymin
							|| triangles[input[k]].v1.z == root->box.zmax && triangles[input[k]].v2.z == root->box.zmax && triangles[input[k]].v3.z == root->box.zmax)
						{
							if (triangles[input[k]].v1.x == root->box.xmax && triangles[input[k]].v2.x == root->box.xmax && triangles[input[k]].v3.x == root->box.xmax)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::YZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else if (triangles[input[k]].v1.y == root->box.ymin && triangles[input[k]].v2.y == root->box.ymin && triangles[input[k]].v3.y == root->box.ymin)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XY))
								{
									input_t[i].push_back(input[k]);
								}
							}
						}
						break;
					}
					case 4:
					{
						if (triangles[input[k]].v1.x == xmid && triangles[input[k]].v2.x == xmid && triangles[input[k]].v3.x == xmid)
						{
							tri_exclude[k][5] = 1;
							input_t[i].push_back(input[k]);
						}
						else if (triangles[input[k]].v1.y == ymid && triangles[input[k]].v2.y == ymid && triangles[input[k]].v3.y == ymid)
						{
							tri_exclude[k][7] = 1;
							input_t[i].push_back(input[k]);
						}
						else if (triangles[input[k]].v1.x == root->box.xmax && triangles[input[k]].v2.x == root->box.xmax && triangles[input[k]].v3.x == root->box.xmax
							|| triangles[input[k]].v1.y == root->box.ymax && triangles[input[k]].v2.y == root->box.ymax && triangles[input[k]].v3.y == root->box.ymax
							|| triangles[input[k]].v1.z == root->box.zmin && triangles[input[k]].v2.z == root->box.zmin && triangles[input[k]].v3.z == root->box.zmin)
						{
							if (triangles[input[k]].v1.x == root->box.xmax && triangles[input[k]].v2.x == root->box.xmax && triangles[input[k]].v3.x == root->box.xmax)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::YZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else if (triangles[input[k]].v1.y == root->box.ymax && triangles[input[k]].v2.y == root->box.ymax && triangles[input[k]].v3.y == root->box.ymax)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XY))
								{
									input_t[i].push_back(input[k]);
								}
							}
						}
						break;
					}
					case 5:
					{
						if (triangles[input[k]].v1.y == ymid && triangles[input[k]].v2.y == ymid && triangles[input[k]].v3.y == ymid)
						{
							tri_exclude[k][6] = 1;
							input_t[i].push_back(input[k]);
						}
						else if (triangles[input[k]].v1.x == root->box.xmin && triangles[input[k]].v2.x == root->box.xmin && triangles[input[k]].v3.x == root->box.xmin
							|| triangles[input[k]].v1.y == root->box.ymax && triangles[input[k]].v2.y == root->box.ymax && triangles[input[k]].v3.y == root->box.ymax
							|| triangles[input[k]].v1.z == root->box.zmin && triangles[input[k]].v2.z == root->box.zmin && triangles[input[k]].v3.z == root->box.zmin)
						{
							if (triangles[input[k]].v1.x == root->box.xmin && triangles[input[k]].v2.x == root->box.xmin && triangles[input[k]].v3.x == root->box.xmin)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::YZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else if (triangles[input[k]].v1.y == root->box.ymax && triangles[input[k]].v2.y == root->box.ymax && triangles[input[k]].v3.y == root->box.ymax)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XY))
								{
									input_t[i].push_back(input[k]);
								}
							}
						}
						break;
					}
					case 6:
					{
						if (triangles[input[k]].v1.x == xmid && triangles[input[k]].v2.x == xmid && triangles[input[k]].v3.x == xmid)
						{
							tri_exclude[k][7] = 1;
							input_t[i].push_back(input[k]);
						}
						else if (triangles[input[k]].v1.x == root->box.xmin && triangles[input[k]].v2.x == root->box.xmin && triangles[input[k]].v3.x == root->box.xmin
							|| triangles[input[k]].v1.y == root->box.ymin && triangles[input[k]].v2.y == root->box.ymin && triangles[input[k]].v3.y == root->box.ymin
							|| triangles[input[k]].v1.z == root->box.zmin && triangles[input[k]].v2.z == root->box.zmin && triangles[input[k]].v3.z == root->box.zmin)
						{
							if (triangles[input[k]].v1.x == root->box.xmin && triangles[input[k]].v2.x == root->box.xmin && triangles[input[k]].v3.x == root->box.xmin)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::YZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else if (triangles[input[k]].v1.y == root->box.ymin && triangles[input[k]].v2.y == root->box.ymin && triangles[input[k]].v3.y == root->box.ymin)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XY))
								{
									input_t[i].push_back(input[k]);
								}
							}
						}
						break;
					}
					case 7:
					{
						if (triangles[input[k]].v1.x == root->box.xmax && triangles[input[k]].v2.x == root->box.xmax && triangles[input[k]].v3.x == root->box.xmax
							|| triangles[input[k]].v1.y == root->box.ymin && triangles[input[k]].v2.y == root->box.ymin && triangles[input[k]].v3.y == root->box.ymin
							|| triangles[input[k]].v1.z == root->box.zmin && triangles[input[k]].v2.z == root->box.zmin && triangles[input[k]].v3.z == root->box.zmin)
						{
							if (triangles[input[k]].v1.x == root->box.xmax && triangles[input[k]].v2.x == root->box.xmax && triangles[input[k]].v3.x == root->box.xmax)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::YZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else if (triangles[input[k]].v1.y == root->box.ymin && triangles[input[k]].v2.y == root->box.ymin && triangles[input[k]].v3.y == root->box.ymin)
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XZ))
								{
									input_t[i].push_back(input[k]);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(sub_box[i], triangles[input[k]], InSurface::XY))
								{
									input_t[i].push_back(input[k]);
								}
							}
						}
						break;
					}
					}
				}

				break;
			}
			}	
		}
		input_t[i].shrink_to_fit();
	}

	root->triangles.clear();
	for (size_t i = 0; i < 8; ++i)
	{
		if (input_t[i].empty())
			continue;

		root->child[i] = new OctTreeNode();
		root->child[i]->box = sub_box[i];
		buildOctTree(input_t[i], root->child[i], depth + 1, false);
	}
}

bool OctTree::insert(const Triangle& t)
{
	if (root == nullptr)
	{
		root = new OctTreeNode();
		root->is_leaf = true;
		AABB b;
		b.xmin = min({ t.v1.x, t.v2.x, t.v3.x });
		b.ymin = min({ t.v1.y, t.v2.y, t.v3.y });
		b.zmin = min({ t.v1.z, t.v2.z, t.v3.z });
		b.xmax = max({ t.v1.x, t.v2.x, t.v3.x });
		b.ymax = max({ t.v1.y, t.v2.y, t.v3.y });
		b.zmax = max({ t.v1.z, t.v2.z, t.v3.z });
		root->box = b;
		root->triangles.push_back(0);
		triangles.push_back(t);
		is_removable.push_back(0);
		root->num_of_tri_not_deleted = 1;
		return true;
	}

	size_t i = 0;
	for (; i < triangles.size(); ++i)
	{
		if (triangles[i] == t)
		{
			if (is_removable[i] == 0)
			{
				return false;
			}
			is_removable[i] = 0;
			break;
		}
	}

	bool insert_be_deleted;
	if (i == triangles.size())
	{
		bool r1 = t.v1.x >= root->box.xmin && t.v1.x <= root->box.xmax &&
			t.v1.y >= root->box.ymin && t.v1.y <= root->box.ymax
			&& t.v1.z >= root->box.zmin && t.v1.z <= root->box.zmax;
		bool r2 = t.v2.x >= root->box.xmin && t.v2.x <= root->box.xmax
			&& t.v2.y >= root->box.ymin && t.v2.y <= root->box.ymax
			&& t.v2.z >= root->box.zmin && t.v2.z <= root->box.zmax;
		bool r3 = t.v3.x >= root->box.xmin && t.v3.x <= root->box.xmax
			&& t.v3.y >= root->box.ymin && t.v3.y <= root->box.ymax
			&& t.v3.z >= root->box.zmin && t.v3.z <= root->box.zmax;

		if (!r1 || !r2 || !r3)
			return false;

		triangles.push_back(t);
		is_removable.push_back(0);
		insert_be_deleted = false;
	}
	else
	{
		insert_be_deleted = true;
	}

	insertTriangleToOctTree(root, t, insert_be_deleted, i);
	return true;
}


void OctTree::insertTriangleToOctTree(OctTreeNode* cur, const Triangle& t, bool insert_be_deleted, size_t index)
{
	if (cur->is_leaf == false)
	{
		RationalNumber xmid = (cur->box.xmin + cur->box.xmax) / RationalNumber(2, 1);
		RationalNumber ymid = (cur->box.ymin + cur->box.ymax) / RationalNumber(2, 1);
		RationalNumber zmid = (cur->box.zmin + cur->box.zmax) / RationalNumber(2, 1);
		vector<int> tri_exclude(8, 0);

		for (size_t i = 0; i < 8; ++i)
		{
			if (tri_exclude[i])
				continue;

			AABB b;
			if (cur->child[i] != nullptr)
			{
				b = cur->child[i]->box;
			}
			else
			{
				computeAABB(b, xmid, ymid, zmid, cur->box, i);

			}

			switch (triangleAndAABBIntersectPrecise(b, t))
			{
			case 1:
			{
				break; 
			}
			case 2:
			{
				for (size_t m = i + 1; m < 8; ++m)
				{
					tri_exclude[m] = 1;
				}
				insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
				break;
			}
			case 3:   
			{
				if (triangleHasPartStrictInAABB(b, t))  
				{
					insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
				}
				else
				{
					switch (i)  
					{  
					case 0: 
					{
						if (t.v1.x == xmid && t.v2.x == xmid && t.v3.x == xmid)
						{
							tri_exclude[1] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
							continue;
						}
						else if (t.v1.y == ymid && t.v2.y == ymid && t.v3.y == ymid)
						{
							tri_exclude[3] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
							continue;
						}
						else if (t.v1.z == zmid && t.v2.z == zmid && t.v3.z == zmid)
						{
							tri_exclude[4] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
							continue;
						}
						else if (t.v1.x == cur->box.xmax && t.v2.x == cur->box.xmax && t.v3.x == cur->box.xmax
							|| t.v1.y == cur->box.ymax && t.v2.y == cur->box.ymax && t.v3.y == cur->box.ymax
							|| t.v1.z == cur->box.zmax && t.v2.z == cur->box.zmax && t.v3.z == cur->box.zmax)
						{
							if (t.v1.x == cur->box.xmax && t.v2.x == cur->box.xmax && t.v3.x == cur->box.xmax)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::YZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else if (t.v1.y == cur->box.ymax && t.v2.y == cur->box.ymax && t.v3.y == cur->box.ymax)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XY))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
						}
						break;
					}
					case 1:
					{
						if (t.v1.y == ymid && t.v2.y == ymid && t.v3.y == ymid)
						{
							tri_exclude[2] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
						}
						else if (t.v1.z == zmid && t.v2.z == zmid && t.v3.z == zmid)
						{
							tri_exclude[5] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
						}
						else if (t.v1.x == cur->box.xmin && t.v2.x == cur->box.xmin && t.v3.x == cur->box.xmin
							|| t.v1.y == cur->box.ymax && t.v2.y == cur->box.ymax && t.v3.y == cur->box.ymax
							|| t.v1.z == cur->box.zmax && t.v2.z == cur->box.zmax && t.v3.z == cur->box.zmax)
						{
							if (t.v1.x == cur->box.xmin && t.v2.x == cur->box.xmin && t.v3.x == cur->box.xmin)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::YZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else if (t.v1.y == cur->box.ymax && t.v2.y == cur->box.ymax && t.v3.y == cur->box.ymax)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XY))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
						}
						break;
					}
					case 2:
					{
						if (t.v1.x == xmid && t.v2.x == xmid && t.v3.x == xmid)
						{
							tri_exclude[3] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
						}
						else if (t.v1.z == zmid && t.v2.z == zmid && t.v3.z == zmid)
						{
							tri_exclude[6] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
						}
						else if (t.v1.x == cur->box.xmin && t.v2.x == cur->box.xmin && t.v3.x == cur->box.xmin
							|| t.v1.y == cur->box.ymin && t.v2.y == cur->box.ymin && t.v3.y == cur->box.ymin
							|| t.v1.z == cur->box.zmax && t.v2.z == cur->box.zmax && t.v3.z == cur->box.zmax)
						{
							if (t.v1.x == cur->box.xmin && t.v2.x == cur->box.xmin && t.v3.x == cur->box.xmin)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::YZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else if (t.v1.y == cur->box.ymin && t.v2.y == cur->box.ymin && t.v3.y == cur->box.ymin)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XY))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
						}
						break;
					}
					case 3:
					{
						if (t.v1.z == zmid && t.v2.z == zmid && t.v3.z == zmid)
						{
							tri_exclude[7] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
						}
						else if (t.v1.x == cur->box.xmax && t.v2.x == cur->box.xmax && t.v3.x == cur->box.xmax
							|| t.v1.y == cur->box.ymin && t.v2.y == cur->box.ymin && t.v3.y == cur->box.ymin
							|| t.v1.z == cur->box.zmax && t.v2.z == cur->box.zmax && t.v3.z == cur->box.zmax)
						{
							if (t.v1.x == cur->box.xmax && t.v2.x == cur->box.xmax && t.v3.x == cur->box.xmax)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::YZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else if (t.v1.y == cur->box.ymin && t.v2.y == cur->box.ymin && t.v3.y == cur->box.ymin)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XY))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
						}
						break;
					}
					case 4:
					{
						if (t.v1.x == xmid && t.v2.x == xmid && t.v3.x == xmid)
						{
							tri_exclude[5] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
						}
						else if (t.v1.y == ymid && t.v2.y == ymid && t.v3.y == ymid)
						{
							tri_exclude[7] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
						}
						else if (t.v1.x == cur->box.xmax && t.v2.x == cur->box.xmax && t.v3.x == cur->box.xmax
							|| t.v1.y == cur->box.ymax && t.v2.y == cur->box.ymax && t.v3.y == cur->box.ymax
							|| t.v1.z == cur->box.zmin && t.v2.z == cur->box.zmin && t.v3.z == cur->box.zmin)
						{
							if (t.v1.x == cur->box.xmax && t.v2.x == cur->box.xmax && t.v3.x == cur->box.xmax)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::YZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else if (t.v1.y == cur->box.ymax && t.v2.y == cur->box.ymax && t.v3.y == cur->box.ymax)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XY))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
						}
						break;
					}
					case 5:
					{
						if (t.v1.y == ymid && t.v2.y == ymid && t.v3.y == ymid)
						{
							tri_exclude[6] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
						}
						else if (t.v1.x == cur->box.xmin && t.v2.x == cur->box.xmin && t.v3.x == cur->box.xmin
							|| t.v1.y == cur->box.ymax && t.v2.y == cur->box.ymax && t.v3.y == cur->box.ymax
							|| t.v1.z == cur->box.zmin && t.v2.z == cur->box.zmin && t.v3.z == cur->box.zmin)
						{
							if (t.v1.x == cur->box.xmin && t.v2.x == cur->box.xmin && t.v3.x == cur->box.xmin)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::YZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else if (t.v1.y == cur->box.ymax && t.v2.y == cur->box.ymax && t.v3.y == cur->box.ymax)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XY))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
						}
						break;
					}
					case 6:
					{
						if (t.v1.x == xmid && t.v2.x == xmid && t.v3.x == xmid)
						{
							tri_exclude[7] = 1;
							insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
						}
						else if (t.v1.x == cur->box.xmin && t.v2.x == cur->box.xmin && t.v3.x == cur->box.xmin
							|| t.v1.y == cur->box.ymin && t.v2.y == cur->box.ymin && t.v3.y == cur->box.ymin
							|| t.v1.z == cur->box.zmin && t.v2.z == cur->box.zmin && t.v3.z == cur->box.zmin)
						{
							if (t.v1.x == cur->box.xmin && t.v2.x == cur->box.xmin && t.v3.x == cur->box.xmin)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::YZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else if (t.v1.y == cur->box.ymin && t.v2.y == cur->box.ymin && t.v3.y == cur->box.ymin)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XY))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
						}
						break;
					}
					case 7:
					{
						if (t.v1.x == cur->box.xmax && t.v2.x == cur->box.xmax && t.v3.x == cur->box.xmax
							|| t.v1.y == cur->box.ymin && t.v2.y == cur->box.ymin && t.v3.y == cur->box.ymin
							|| t.v1.z == cur->box.zmin && t.v2.z == cur->box.zmin && t.v3.z == cur->box.zmin)
						{
							if (t.v1.x == cur->box.xmax && t.v2.x == cur->box.xmax && t.v3.x == cur->box.xmax)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::YZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else if (t.v1.y == cur->box.ymin && t.v2.y == cur->box.ymin && t.v3.y == cur->box.ymin)
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XZ))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
							else
							{
								if (triangleInterAABBWithOneSurface(b, t, InSurface::XY))
								{
									insertIntoSubTree(cur, i, b, index, t, insert_be_deleted);
								}
							}
						}
						break;
					}
					}
				}

				break;
			}
			}
		}
		return;
	}

	++(cur->num_of_tri_not_deleted);
	if (insert_be_deleted == false)
	{
		cur->triangles.push_back(index);
	}
	else
	{
		auto r = lower_bound(cur->triangles.begin(), cur->triangles.end(), index);
		if (r == cur->triangles.end() || index != *r)
		{
			cur->triangles.insert(r, index);
		}
	}

	buildOctTree(cur->triangles, cur, 1, true);
}

void OctTree::insertIntoSubTree(OctTreeNode* cur, size_t i, const AABB& b, size_t index, const Triangle& t, bool insert_be_deleted)
{
	if (cur->child[i] == nullptr)
	{
		cur->child[i] = new OctTreeNode();
		cur->child[i]->is_leaf = true;
		cur->child[i]->box = b;
		cur->child[i]->triangles.push_back(index);
		cur->child[i]->num_of_tri_not_deleted = 1;
	}
	else
	{
		insertTriangleToOctTree(cur->child[i], t, insert_be_deleted, index);
	}
}

int main()
{
	using std::vector;
	Triangle t1;
	Triangle t2;
	Triangle t3;
	t1.v1.x = RationalNumber(0, 1);
	t1.v1.y = RationalNumber(0, 1);
	t1.v1.z = RationalNumber(0, 1);
	t1.v2.x = RationalNumber(1, 1);
	t1.v2.y = RationalNumber(0, 1);
	t1.v2.z = RationalNumber(0, 1);
	t1.v3.x = RationalNumber(0, 1);
	t1.v3.y = RationalNumber(1, 1);
	t1.v3.z = RationalNumber(0, 1);

	t2.v1.x = RationalNumber(0, 1);
	t2.v1.y = RationalNumber(0, 1);
	t2.v1.z = RationalNumber(1, 1);
	t2.v2.x = RationalNumber(1, 1);
	t2.v2.y = RationalNumber(0, 1);
	t2.v2.z = RationalNumber(1, 1);
	t2.v3.x = RationalNumber(0, 1);
	t2.v3.y = RationalNumber(1, 1);
	t2.v3.z = RationalNumber(1, 1);

	t3.v1.x = RationalNumber(2, 1);
	t3.v1.y = RationalNumber(0, 1);
	t3.v1.z = RationalNumber(0, 1);
	t3.v2.x = RationalNumber(3, 1);
	t3.v2.y = RationalNumber(0, 1);
	t3.v2.z = RationalNumber(0, 1);
	t3.v3.x = RationalNumber(2, 1);
	t3.v3.y = RationalNumber(1, 1);
	t3.v3.z = RationalNumber(0, 1);

	vector<Triangle> test_tri(3);
	test_tri[0] = t1;
	test_tri[1] = t2;
	test_tri[2] = t3;

	OctTree test_tree(test_tri, 5, 1);
	test_tree.removeTriangle(2);
	test_tree.insert(t3);

	PointOrVector3D _1p;
	_1p.x = RationalNumber(1, 2);
	_1p.y = RationalNumber(1, 2);
	_1p.z = RationalNumber(-1, 1);
	PointOrVector3D _1d;
	_1d.x = RationalNumber(0, 1);
	_1d.y = RationalNumber(0, 1);
	_1d.z = RationalNumber(1, 1);

	pair<std::vector<pair<size_t, InterType>>, Value> r1 = test_tree.solveLineAndTriangleInter(_1p, _1d, Value(r_zero), Value(v_type::positive_inf));
	assert(r1.second != v_type::positive_inf && r1.first.size() == 1 && test_tree.getTriangles()[r1.first[0].first] == t1 && r1.second.type == v_type::finite && r1.second.value == RationalNumber(1, 1) && "_1faild!");

	PointOrVector3D _2p;
	_2p.x = RationalNumber(1, 2);
	_2p.y = RationalNumber(1, 2);
	_2p.z = RationalNumber(1, 2);
	PointOrVector3D _2d;
	_2d.x = RationalNumber(0, 1);
	_2d.y = RationalNumber(0, 1);
	_2d.z = RationalNumber(1, 1);

	pair<std::vector<pair<size_t, InterType>>, Value> r2 = test_tree.solveLineAndTriangleInter(_2p, _2d, Value(r_zero), Value(v_type::positive_inf));
	assert(r2.second != v_type::positive_inf && r2.first.size() == 1 && test_tree.getTriangles()[r2.first[0].first] == t2 && r2.second.type == v_type::finite && r2.second.value == RationalNumber(1, 2) && "_2faild!");

	PointOrVector3D _3p;
	_3p.x = RationalNumber(5, 2);
	_3p.y = RationalNumber(1, 2);
	_3p.z = RationalNumber(-1, 1);
	PointOrVector3D _3d;
	_3d.x = RationalNumber(0, 1);
	_3d.y = RationalNumber(0, 1);
	_3d.z = RationalNumber(1, 1);

	pair<std::vector<pair<size_t, InterType>>, Value> r3 = test_tree.solveLineAndTriangleInter(_3p, _3d, Value(r_zero), Value(v_type::positive_inf));
	assert(r3.second != v_type::positive_inf && r3.first.size() == 1 && test_tree.getTriangles()[r3.first[0].first] == t3 && r3.second.type == v_type::finite && r3.second.value == RationalNumber(1, 1) && "_3faild!");

	PointOrVector3D _4p;
	_4p.x = RationalNumber(1, 2);
	_4p.y = RationalNumber(-1, 1);
	_4p.z = RationalNumber(0, 1);
	PointOrVector3D _4d;
	_4d.x = RationalNumber(0, 1);
	_4d.y = RationalNumber(1, 1);
	_4d.z = RationalNumber(0, 1);

	pair<std::vector<pair<size_t, InterType>>, Value> r4 = test_tree.solveLineAndTriangleInter(_4p, _4d, Value(r_zero), Value(v_type::positive_inf));
	assert(r4.second == v_type::positive_inf && "_4faild!");

	PointOrVector3D _5p;
	_5p.x = RationalNumber(1, 2);
	_5p.y = RationalNumber(1, 2);
	_5p.z = RationalNumber(2, 1);
	PointOrVector3D _5d;
	_5d.x = RationalNumber(0, 1);
	_5d.y = RationalNumber(0, 1);
	_5d.z = RationalNumber(-1, 1);

	pair<std::vector<pair<size_t, InterType>>, Value> r5 = test_tree.solveLineAndTriangleInter(_5p, _5d, Value(r_zero), Value(v_type::positive_inf));
	assert(r5.second != v_type::positive_inf && r5.first.size() == 1 && test_tree.getTriangles()[r5.first[0].first]== t2 && r5.second.type == v_type::finite && r5.second.value == RationalNumber(1, 1) && "_5faild!");

	PointOrVector3D _6p;
	_6p.x = RationalNumber(10, 1);
	_6p.y = RationalNumber(10, 1);
	_6p.z = RationalNumber(10, 1);
	PointOrVector3D _6d;
	_6d.x = RationalNumber(-1, 1);
	_6d.y = RationalNumber(-1, 1);
	_6d.z = RationalNumber(-1, 1);

	pair<std::vector<pair<size_t, InterType>>, Value> r6 = test_tree.solveLineAndTriangleInter(_6p, _6d, Value(r_zero), Value(v_type::positive_inf));
	assert(r6.second != v_type::positive_inf && r6.first.size() == 1 && test_tree.getTriangles()[r6.first[0].first] == t1 && r6.second.type == v_type::finite && r6.second.value == RationalNumber(10, 1) && "_6faild!");

	PointOrVector3D _7p;
	_7p.x = RationalNumber(1, 2);
	_7p.y = RationalNumber(1, 2);
	_7p.z = RationalNumber(1, 2);
	PointOrVector3D _7d;
	_7d.x = RationalNumber(1, 1);
	_7d.y = RationalNumber(0, 1);
	_7d.z = RationalNumber(0, 1);

	pair<std::vector<pair<size_t, InterType>>, Value> r7 = test_tree.solveLineAndTriangleInter(_7p, _7d, Value(r_zero), Value(v_type::positive_inf));
	assert(r7.second == v_type::positive_inf && "_7faild!");

	PointOrVector3D _8p;
	_8p.x = RationalNumber(4, 5);
	_8p.y = RationalNumber(4, 5);
	_8p.z = RationalNumber(-1, 1);
	PointOrVector3D _8d;
	_8d.x = RationalNumber(0, 1);
	_8d.y = RationalNumber(0, 1);
	_8d.z = RationalNumber(1, 1);

	pair<std::vector<pair<size_t, InterType>>, Value> r8 = test_tree.solveLineAndTriangleInter(_8p, _8d, Value(r_zero), Value(v_type::positive_inf));
	assert(r8.second == v_type::positive_inf && "_8faild!");

	PointOrVector3D _9p;
	_9p.x = RationalNumber(4, 5);
	_9p.y = RationalNumber(4, 5);
	_9p.z = RationalNumber(1, 2);
	PointOrVector3D _9d;
	_9d.x = RationalNumber(0, 1);
	_9d.y = RationalNumber(0, 1);
	_9d.z = RationalNumber(1, 1);

	pair<std::vector<pair<size_t, InterType>>, Value> r9 = test_tree.solveLineAndTriangleInter(_9p, _9d, Value(r_zero), Value(v_type::positive_inf));
	assert(r9.second == v_type::positive_inf && "_9faild!");
    return 0;
}

