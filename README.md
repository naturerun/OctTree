# 概述
这是一个利用八叉树实现光线追踪中光线和三角形加速求交的库
RationalNumber.h是有理数类的实现
utility.h是主程序需要用到的工具函数
八叉树.cpp是主程序文件
八叉树实现能够分割空间，高效存储和索引三角形,并且允许三角形的动态插入和删除
使用示例:
```cpp
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
```
t1,t2,t3为Triangle类型的对象，代表一个三角形,Triangle类型的成员即为构成三角形的三个顶点,这里不再赘述
带有三个参数的OctTree构造函数的第一个参数是保存在即将创建的八叉树中的三角形列表,第二个参数是空间分割的最大递归深度,第三个参数是三角形容量,限制八叉树叶节点的三角形数量上限，超出该上限值叶节点需要分割
带有两个参数的构造函数创建初始化一棵空的八叉树
八叉树类的成员triangles为保存在八叉树中的所有三角形集合,八叉树用triangles的索引保存三角形。is_removable是删除标记，指明triangles中的三角形是否已经删除
removeTriangle成员函数删除triangles的索引为index的三角形
insert成员函数向八叉树插入一个三角形
solveLineAndTriangleInter返回和起始点在p0方向向量为d,t在[tmin,tmax]区间的线段相交的所有三角形
具体用法参看主程序文件的main函数