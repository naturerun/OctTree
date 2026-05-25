#pragma once
using std::pair;

RationalNumber r_zero(0, 1);

enum class v_type { positive_inf, negative_inf, finite };

struct Value
{
	RationalNumber value;
	v_type type;
	Value() = default;
	Value(RationalNumber _v) :value(_v), type(v_type::finite) {}
	Value(v_type v_t) :type(v_t), value(0, 1) {}
	bool operator<(const Value& be) const
	{
		if (type == be.type)
		{
			if (type == v_type::finite)
			{
				return value < be.value;
			}
			return false;
		}
		return type == v_type::negative_inf && (be.type == v_type::finite || be.type == v_type::positive_inf)
			|| type == v_type::finite && be.type == v_type::positive_inf;
	}

	Value operator-()
	{
		if (type == v_type::negative_inf)
		{
			return Value(v_type::positive_inf);
		}

		if (type == v_type::positive_inf)
		{
			return Value(v_type::negative_inf);
		}

		return Value(-value);
	}

	Value operator-(Value& be)
	{
		return (*this) + (-be);
	}

	Value operator+(const Value& be)
	{
		Value zero{ r_zero };
		if (type == v_type::negative_inf)
		{
			if (be.type == v_type::positive_inf)
			{
				return zero;
			}
			return Value(v_type::negative_inf);
		}

		else if (type == v_type::positive_inf)
		{
			if (be.type == v_type::negative_inf)
			{
				return zero;
			}
			return Value(v_type::positive_inf);
		}
		else
		{
			if (be.type == v_type::negative_inf)
			{
				return Value(v_type::negative_inf);
			}
			else if (be.type == v_type::positive_inf)
			{
				return Value(v_type::positive_inf);
			}
			else
			{
				return Value(value + be.value);
			}
		}
	}

	Value operator*(const Value& be)
	{
		const Value zero{ r_zero };
		if (be == zero || (*this) == zero)
			return zero;

		if (be.type != type)
		{
			if (type == v_type::negative_inf)
			{
				if (be.type == v_type::finite)
				{
					if (be > zero)
					{
						return Value(v_type::negative_inf);

					}
					return Value(v_type::positive_inf);
				}
				return Value(v_type::negative_inf);
			}

			if (type == v_type::positive_inf)
			{
				if (be.type == v_type::finite)
				{
					if (be > zero)
					{
						return Value(v_type::positive_inf);

					}
					return Value(v_type::negative_inf);
				}
				return Value(v_type::negative_inf);
			}

			if ((*this) > zero)
			{
				if (be.type == v_type::positive_inf)
				{
					return Value(v_type::positive_inf);
				}
				return Value(v_type::negative_inf);
			}
			else
			{
				if (be.type == v_type::positive_inf)
				{
					return Value(v_type::negative_inf);
				}
				return Value(v_type::positive_inf);
			}
		}

		if (type == v_type::finite)
		{
			return Value(value * be.value);
		}
		return Value(v_type::positive_inf);
	}

	bool operator==(const Value& be) const
	{
		if (type == be.type)
		{
			if (type == v_type::finite)
			{
				return value == be.value;
			}
			return true;
		}
		return false;
	}

	bool operator!=(const Value& be) const
	{
		return !operator == (be);
	}

	bool operator<=(const Value& be) const
	{
		return !operator>(be);
	}

	bool operator>=(const Value& be) const
	{
		return !operator<(be);
	}

	bool operator>(const Value& be) const
	{
		return operator <(be) == false && operator==(be) == false;
	}
};

Value zero{ r_zero };
pair<Value, Value> getIntersection(pair<Value, Value> l, pair<Value, Value> j, bool& inter) //返回区间l,j的交区间,l,g左右端点均为正整数,返回(0,0)表示交集不存在
{
	if (l.second < j.first || l.first > j.second)
	{
		inter = false;
		return { zero , zero };
	}

	inter = true;
	if (j.first <= l.first)
	{
		if (j.second < l.second)
			return { l.first, j.second };
		return l;
	}

	if (j.second <= l.second)
		return j;
	return { j.first, l.second };
}

RationalNumber compute_t(RationalNumber c, RationalNumber C, RationalNumber value, bool& inf)
{
	if (C == r_zero)
	{
		inf = true;
		return r_zero;
	}

	inf = false;
	return (value - c) / C;
}
//(xmin, ymin, zmin)是AABB下侧矩形左下角的顶点,(xmax, ymax, zmax)是AABB上侧矩形右上角的顶点,线段为(x0,y0,z0)+t*(A,B,C),t在tmin,tmax之间(包括tmin,tmax),tmin可以等于小于大于tmax
//当返回值的first成员为v_type::positive_inf时,线段和AABB不相交,否则为线段和AABB相交部分的参数t变化范围[first,second],返回的pair满足first<=second,first和second均为有限值
pair<Value, Value> judgeAABBInterline(RationalNumber xmin, RationalNumber ymin, RationalNumber zmin, RationalNumber xmax, RationalNumber ymax, RationalNumber zmax,
	RationalNumber x0, RationalNumber y0, RationalNumber z0, RationalNumber A, RationalNumber B, RationalNumber C, RationalNumber tmin, RationalNumber tmax,
	v_type tmin_type, v_type tmax_type)
{
	Value _x0{ x0 };
	Value _y0{ y0 };
	Value _z0{ z0 };
	Value _A{ A };
	Value _B{ B };
	Value _C{ C };
	Value _xmin{ xmin };
	Value _xmax{ xmax };
	Value _ymin{ ymin };
	Value _ymax{ ymax };
	Value _zmin{ zmin };
	Value _zmax{ zmax };

	Value _tmin;
	if (tmin_type == v_type::finite)
	{
		_tmin = Value{ tmin };
	}
	else
	{
		_tmin = Value{ tmin_type };
	}

	Value _tmax;
	if (tmax_type == v_type::finite)
	{
		_tmax = Value{ tmax };
	}
	else
	{
		_tmax = Value{ tmax_type };
	}

	Value x_range_left = (_x0 + _A * _tmin);
	Value x_range_right = (_x0 + _A * _tmax);
	Value y_range_left = (_y0 + _B * _tmin);
	Value y_range_right = (_y0 + _B * _tmax);
	Value z_range_left = (_z0 + _C * _tmin);
	Value z_range_right = (_z0 + _C * _tmax);

	bool x_exchange = false;
	if (x_range_left > x_range_right)
	{
		std::swap(x_range_left, x_range_right);
		x_exchange = true;
	}

	bool y_exchange = false;
	if (y_range_left > y_range_right)
	{
		std::swap(y_range_left, y_range_right);
		y_exchange = true;
	}

	bool z_exchange = false;
	if (z_range_left > z_range_right)
	{
		std::swap(z_range_left, z_range_right);
		z_exchange = true;
	}

	bool x_is_inter;
	pair<Value, Value> x_inter = getIntersection({ x_range_left, x_range_right }, { _xmin, _xmax }, x_is_inter);
	if (x_is_inter == false)
		return { Value{v_type::positive_inf}, Value{v_type::positive_inf} };

	bool y_is_inter;
	pair<Value, Value> y_inter = getIntersection({ y_range_left, y_range_right }, { _ymin, _ymax }, y_is_inter);
	if (y_is_inter == false)
		return { Value{v_type::positive_inf}, Value{v_type::positive_inf} };

	bool z_is_inter;
	pair<Value, Value> z_inter = getIntersection({ z_range_left, z_range_right }, { _zmin, _zmax }, z_is_inter);
	if (z_is_inter == false)
		return { Value{v_type::positive_inf}, Value{v_type::positive_inf} };

	Value t_x_min;
	if (x_inter.first != x_range_left)
	{
		bool inf;
		t_x_min = Value(compute_t(_x0.value, _A.value, x_inter.first.value, inf));
	}
	else
	{
		t_x_min = _tmin;
	}

	Value t_x_max;
	if (x_inter.second != x_range_right)
	{
		bool inf;
		t_x_max = Value(compute_t(_x0.value, _A.value, x_inter.second.value, inf));
	}
	else
	{
		t_x_max = _tmax;
	}

	if (t_x_min > t_x_max)
	{
		std::swap(t_x_min, t_x_max);
	}


	Value t_y_min;
	if (y_inter.first != y_range_left)
	{
		bool inf;
		t_y_min = Value(compute_t(_y0.value, _B.value, y_inter.first.value, inf));
	}
	else
	{
		t_y_min = _tmin;
	}

	Value t_y_max;
	if (y_inter.second != y_range_right)
	{
		bool inf;
		t_y_max = Value(compute_t(_y0.value, _B.value, y_inter.second.value, inf));
	}
	else
	{
		t_y_max = _tmax;
	}

	if (t_y_min > t_y_max)
	{
		std::swap(t_y_min, t_y_max);
	}

	Value t_z_min;
	if (z_inter.first != z_range_left)
	{
		bool inf;
		t_z_min = Value(compute_t(_z0.value, _C.value, z_inter.first.value, inf));
	}
	else
	{
		t_z_min = _tmin;
	}

	Value t_z_max;
	if (z_inter.second != z_range_right)
	{
		bool inf;
		t_z_max = Value(compute_t(_z0.value, _C.value, z_inter.second.value, inf));
	}
	else
	{
		t_z_max = _tmax;
	}

	if (t_z_min > t_z_max)
	{
		std::swap(t_z_min, t_z_max);
	}

	bool inter_xy;
	pair<Value, Value> xy = getIntersection({ t_x_min, t_x_max }, { t_y_min, t_y_max }, inter_xy);
	if (inter_xy == false)
		return { Value{v_type::positive_inf}, Value{v_type::positive_inf} };

	bool inter_xyz;
	pair<Value, Value> xyz = getIntersection(xy, { t_z_min, t_z_max }, inter_xyz);
	if (inter_xyz == false)
		return { Value{v_type::positive_inf}, Value{v_type::positive_inf} };
	return xyz;
}

struct PointOrVector3D
{
	RationalNumber x;
	RationalNumber y;
	RationalNumber z;
	PointOrVector3D() = default;
	PointOrVector3D(RationalNumber x, RationalNumber y, RationalNumber z) :x(x), y(y), z(z) {}
	bool operator==(const PointOrVector3D& be) const
	{
		return x == be.x && y == be.y && z == be.z;
	}
	PointOrVector3D operator-(const PointOrVector3D &op) const
	{
		PointOrVector3D res;
		res.x = x - op.x;
		res.y = y - op.y;
		res.z = z - op.z;
		return res;
	}
};

RationalNumber innerProduct(PointOrVector3D& v1, PointOrVector3D& v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

PointOrVector3D crossProduct(PointOrVector3D& v1, PointOrVector3D& v2)
{
	PointOrVector3D res;
	res.x = v1.y * v2.z - v1.z * v2.y;
	res.y = v1.z * v2.x - v1.x * v2.z;
	res.z = v1.x * v2.y - v1.y * v2.x;
	return res;
}


PointOrVector3D subtract(PointOrVector3D& v1, PointOrVector3D& v2)
{
	PointOrVector3D res;
	res.x = v1.x - v2.x;
	res.y = v1.y - v2.y;
	res.z = v1.z - v2.z;
	return res;
}

enum class InterType { V1, V3, V2, V1V3, V3V2, V2V1, Internal, NOInter }; //三角形从V1出发逆时针依次为v1,v3,v2
pair<RationalNumber, InterType> lineInterTriangle(PointOrVector3D& p0, PointOrVector3D& d, Value tmin, Value tmax, PointOrVector3D& v1, PointOrVector3D& v2, PointOrVector3D& v3)
{
	PointOrVector3D v1v2 = subtract(v1, v2);
	PointOrVector3D v1v3 = subtract(v1, v3);
	PointOrVector3D n = crossProduct(v1v3, v1v2);
	if (innerProduct(n, d) == r_zero)  //线段在平面内或和平面平形均视为光线未击中平面
	{
		return { r_zero, InterType::NOInter };
	}

	PointOrVector3D s = subtract(v1, p0);
	RationalNumber l = innerProduct(n, s) / innerProduct(n, d);
	Value _l(l);

	if (tmin > tmax)
	{
		std::swap(tmin, tmax);
	}

	if (tmin > _l || _l > tmax)
	{
		return { r_zero, InterType::NOInter };
	}

	PointOrVector3D res;
	res.x = p0.x + d.x * l;
	res.y = p0.y + d.y * l;
	res.z = p0.y + d.y * l;

	PointOrVector3D v_v3v1 = subtract(v3, v1);
	PointOrVector3D v_v2v3 = subtract(v2, v3);
	PointOrVector3D v_v1v2 = subtract(v1, v2);
	PointOrVector3D v_rv1 = subtract(res, v1);
	PointOrVector3D v_rv3 = subtract(res, v3);
	PointOrVector3D v_rv2 = subtract(res, v2);
	PointOrVector3D c_1 = crossProduct(v_v3v1, v_rv1);
	PointOrVector3D c_2 = crossProduct(v_v2v3, v_rv3);
	PointOrVector3D c_3 = crossProduct(v_v1v2, v_rv2);

	RationalNumber _c_1;
	RationalNumber _c_2;
	RationalNumber _c_3;
	if ((_c_1 = innerProduct(c_1, n)) >= r_zero && (_c_2 = innerProduct(c_2, n)) >= r_zero && (_c_3 = innerProduct(c_3, n)) >= r_zero)
	{
		if (_c_1 == r_zero)
		{
			if (res.x == v1.x && res.y == v1.y && res.z == v1.z)
			{
				return { l, InterType::V1 };
			}
			else if (res.x == v3.x && res.y == v3.y && res.z == v3.z)
			{
				return { l, InterType::V3 };
			}
			else
			{
				return { l, InterType::V1V3 };
			}
		}

		if (_c_2 == r_zero)
		{
			if (res.x == v2.x && res.y == v2.y && res.z == v2.z)
			{
				return { l, InterType::V2 };
			}
			else
			{
				return { l, InterType::V3V2 };
			}
		}

		if (_c_3 == r_zero)
		{
			return { l, InterType::V2V1 };
		}

		return { l, InterType::Internal };
	}
	return { r_zero, InterType::NOInter };
}

struct Triangle
{
	PointOrVector3D v1;
	PointOrVector3D v2;
	PointOrVector3D v3;
	bool operator==(const Triangle& be) const
	{
		return v1 == be.v1 && v2 == be.v2 && v3 == be.v3;
	}
};


struct AABB
{
	RationalNumber xmin;
	RationalNumber ymin;
	RationalNumber zmin;
	RationalNumber xmax;
	RationalNumber ymax;
	RationalNumber zmax;
};

struct PointOrVector2D
{
	RationalNumber x;
	RationalNumber y;
	PointOrVector2D() = default;
	bool operator==(const PointOrVector2D& be) const
	{
		return x == be.x && y == be.y;
	}
	PointOrVector2D operator-(const PointOrVector2D& op) const
	{
		PointOrVector2D res;
		res.x = x - op.x;
		res.y = y - op.y;
		return res;
	}

	PointOrVector2D(RationalNumber x, RationalNumber y) :x(x), y(y) {}
};

RationalNumber scalarProduct(RationalNumber x1, RationalNumber y1, RationalNumber x2, RationalNumber y2)
{
	return x1 * x2 + y1 * y2;
}

RationalNumber vectorProduct(RationalNumber x1, RationalNumber y1, RationalNumber x2, RationalNumber y2)
{
	return x1 * y2 - x2 * y1;
}

bool left_turn(RationalNumber x1, RationalNumber y1, RationalNumber x2, RationalNumber y2)
{
	if (vectorProduct(x1, y1, x2, y2) < r_zero)
	{
		return true;
	}
	return false;
}

bool compare(const PointOrVector2D& left_down_point, const PointOrVector2D& left, const PointOrVector2D& right)
{
	if (vectorProduct(left.x - left_down_point.x, left.y - left_down_point.y, right.x - left_down_point.x, right.y - left_down_point.y) > r_zero)
	{
		return true;
	}

	if (vectorProduct(left.x - left_down_point.x, left.y - left_down_point.y, right.x - left_down_point.x, right.y - left_down_point.y) == r_zero)
	{
		if (scalarProduct(left.x - left_down_point.x, left.y - left_down_point.y, left.x - right.x, left.y - right.y) > r_zero)
		{
			return true;
		}
	}
	return false;
}

using std::vector;
using std::swap;
void QuickSort(vector<PointOrVector2D>& seq, int left, int right, const PointOrVector2D& left_down_point)
{
	if (left >= right)
		return;
	int mid = (left + right) / 2;

	if (compare(left_down_point, seq[left], seq[mid]))
		mid = left;
	if (compare(left_down_point, seq[right], seq[mid]))
		mid = right;
	if (mid != left)
	{
		swap(seq[mid], seq[left]);
	}

	if (compare(left_down_point, seq[mid], seq[right]))
	{
		swap(seq[mid], seq[right]);
	}

	int i = left;
	int j = right;
	PointOrVector2D pivot = seq[right];
	while (i < j)
	{
		for (; i < j && compare(left_down_point, pivot, seq[i]) == false; )
		{
			++i;
		}
		seq[j] = seq[i];
		for (; i < j && compare(left_down_point, seq[j], pivot) == false; )
		{
			--j;
		}
		seq[i] = seq[j];
	}
	seq[i] = pivot;
	QuickSort(seq, left, i - 1, left_down_point);
	QuickSort(seq, i + 1, right, left_down_point);
}

void doGrahamScan(vector<PointOrVector2D>& sorted_point, vector<int>& work_stack)
{
	work_stack.push_back(0);
	work_stack.push_back(1);
	for (int i = 2; i < sorted_point.size(); ++i)
	{
		while (left_turn(sorted_point[i].x - sorted_point[work_stack.back()].x, sorted_point[i].y - sorted_point[work_stack.back()].y,
			sorted_point[work_stack.back()].x - sorted_point[work_stack[work_stack.size() - 2]].x, sorted_point[work_stack.back()].y - sorted_point[work_stack[work_stack.size() - 2]].y) == false)
		{
			work_stack.pop_back();
		}
		work_stack.push_back(i);
	}
}

vector<PointOrVector2D> getConvexHull(const vector<PointOrVector2D> &input)  //凸包逆时针方向排列
{
	RationalNumber y_min;
	for (size_t i = 0; i < input.size(); ++i)
	{
		if (i == 0 || input[i].y < y_min)
		{
			y_min = input[i].y;
		}
	}

	RationalNumber x_min;
	int min_x_min_y_index = -1;
	for (size_t i = 0; i < input.size(); ++i)       //这两个循环可以合并为一个循环，可修改
	{
		if (input[i].y == y_min)
		{
			if (min_x_min_y_index == -1 || input[i].x < x_min)
			{
				x_min = input[i].x;
				min_x_min_y_index = i;
			}
		}
	}

	vector<int> work_stack;
	vector<PointOrVector2D> point_fartest;
	vector<PointOrVector2D> other_point;
	other_point.reserve(input.size() - 1);

	for (size_t i = 0; i < input.size(); ++i)
	{
		if (input[i] == input[min_x_min_y_index])
		{
			continue;
		}
		other_point.push_back(input[i]);
	}
	other_point.shrink_to_fit();
	QuickSort(other_point, 0, other_point.size() - 1, { x_min, y_min });
	point_fartest.push_back(other_point[0]);
	for (size_t i = 1; i < other_point.size(); ++i)
	{
		if (vectorProduct(other_point[i - 1].x - x_min, other_point[i - 1].y - y_min, other_point[i].x - x_min, other_point[i].y - y_min) != r_zero)
		{
			point_fartest.push_back(other_point[i]);
		}
	}

	if (point_fartest.size() == 1)
	{
		return { PointOrVector2D(x_min, y_min), PointOrVector2D(point_fartest[0].x, point_fartest[0].y) };
	}
	doGrahamScan(point_fartest, work_stack);
	vector<PointOrVector2D> r_convex_hull(work_stack.size() + 1);
	r_convex_hull[0] = PointOrVector2D(x_min, y_min);
	for (size_t i = 1; i <= work_stack.size(); ++i)
	{
		r_convex_hull[i] = point_fartest[work_stack[i - 1]];
	}
	return r_convex_hull;
}

struct LineSegment
{
	PointOrVector2D left;
	PointOrVector2D right;
	LineSegment(RationalNumber lfirst, RationalNumber lsecond, RationalNumber rfirst, RationalNumber rsecond) :left({ lfirst, lsecond }), right({ rfirst, rsecond }) {}
};

bool pointOnLine(const PointOrVector2D& point, const LineSegment& line)
{
	if (vectorProduct(point.x - line.left.x, point.y - line.left.y, line.right.x - line.left.x, line.right.y - line.left.y) != r_zero)
	{
		return false;
	}

	if (scalarProduct(point.x - line.left.x, point.y - line.left.y, point.x - line.right.x, point.y - line.right.y) > r_zero)
	{
		return false;
	}
	return true;
}

bool parallel(const LineSegment& line1, const LineSegment& line2)
{
	if (vectorProduct(line1.right.x - line1.left.x, line1.right.y - line1.left.y, line2.right.x - line2.left.x, line2.right.y - line2.left.y) == r_zero)
	{
		return true;
	}
	return false;
}

bool cross(const LineSegment& line1, const LineSegment& line2)
{
	RationalNumber vector_product1 = vectorProduct(line1.left.x - line2.left.x, line1.left.y - line2.left.y, line2.right.x - line2.left.x, line2.right.y - line2.left.y);
	RationalNumber vector_product2 = vectorProduct(line1.right.x - line2.left.x, line1.right.y - line2.left.y, line2.right.x - line2.left.x, line2.right.y - line2.left.y);
	if (vector_product1 < r_zero && vector_product2 > r_zero || vector_product1 > r_zero && vector_product2 < r_zero)
	{
		return true;
	}
	return false;
}

bool intersect(const LineSegment& line1, const LineSegment& line2)
{
	if (pointOnLine(line1.left, line2) || pointOnLine(line1.right, line2))
	{
		return true;
	}

	if (pointOnLine(line2.left, line1) || pointOnLine(line2.right, line1))
	{
		return true;
	}

	if (parallel(line1, line2))
	{
		return false;
	}

	if (cross(line1, line2) && cross(line2, line1))
	{
		return true;
	}

	return false;
}