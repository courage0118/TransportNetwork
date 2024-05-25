import csv
import numpy as np
import math
import os


# 1：基本类定义（路段，节点，OD对）
class node():  # 节点类，包含每个节点的坐标信息
    def __init__(self, id=None):
        self.ID = id
        self.origin = -1
        self.innode = []  # 进入节点路段集合
        self.outnode = []  # 离开节点路段集合


class link():  # 路段类，包含路段的起讫点、自由流行驶时间、容量、BRP路阻函数的参数
    def __init__(self):
        self.ID = None
        self.O_link = node()
        self.D_link = node()
        self.FFT = None
        self.Traveltime = None
        self.Capacity = None
        self.B = 0.15  # SiouxFalls网络中路网B与Power系数值均相等，方便期间，统一设置一个初始值
        self.Power = 4


class OD():  # OD类，记录每一个OD需求信息，包含起点、重点、流量需求
    def __init__(self):
        self.ID = None
        self.origin_node = None
        self.Destination = []
        self.odlinks_demand = []


# 2：定义路网核心函数（路网信息导入、UE分配计算等）
class network():  # 网络类，核心函数
    def __init__(self):
        self.Nodes = []  # 网络包含节点集合
        self.Links = []  # 路段集合
        self.Origins = []  # 起点集合，即流量发点
        self.num_Nodes = 0  # 节点个数
        self.num_Links = 0  # 路段个数
        self.num_Origins = 0  # 起点个数
        self.Linkflows = []  # 路段流量集合
        self.Linktimes = []  # 路段行驶时间集合
        self.max_err = 0.001  # UE最大误差
        self.err = 1  # 初始UE误差
        self.Djpathcost = []  # 最短路阻抗集合（Dijkstra算法求得）
        self.Djpath = []  # 最短路集合

    def read_nodes(self, path):
        with open(path, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                newnode = node()
                newnode.ID = int(row[0])
                self.Nodes.append(newnode)

    def read_link(self, path):
        with open(path, 'r') as file:
            reader = csv.reader(file)
            linkid = 1
            for row in reader:
                newlink = link()
                newlink.ID = linkid
                newlink.O_link = node(id=int(row[0]))
                newlink.D_link = node(id=int(row[1]))
                newlink.FFT = int(row[2])
                newlink.Capacity = float(row[3])
                # 将节点信息存储到相应的节点对象中
                newlink.O_link.outnode.append(newlink.ID)
                newlink.D_link.innode.append(newlink.ID)
                linkid += 1
                self.Links.append(newlink)

    def read_od(self, path):
        with open(path, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                newnode = self.Nodes[int(row[0]) - 1]
                if newnode.origin == -1:
                    neworigin = OD()
                    self.num_Origins += 1
                    neworigin.ID = self.num_Origins
                    neworigin.origin_node = newnode
                    self.Nodes[int(row[0]) - 1].origin = neworigin.ID
                    self.Origins.append(neworigin)
                else:
                    neworigin = self.Origins[newnode.origin - 1]
                self.Origins[newnode.origin - 1].Destination.append(self.Nodes[int(row[1]) - 1].ID)
                self.Origins[newnode.origin - 1].odlinks_demand.append(float(row[2]))

    def Dijkstra_path(self, start, end):  # 记录最短路径
        num_nodes = len(self.Nodes)
        unvisited = {node.ID: float('inf') for node in self.Nodes}
        unvisited[start] = 0
        visited = {}
        parents = {node.ID: None for node in self.Nodes}

        while unvisited:
            current_node = min(unvisited, key=unvisited.get)
            current_distance = unvisited[current_node]

            if current_distance == float('inf'):
                break

            for link_id in self.Nodes[current_node - 1].outnode:
                link = self.Links[link_id - 1]
                neighbor = link.D_link.ID
                new_distance = current_distance + link.Traveltime

                if new_distance < unvisited.get(neighbor, float('inf')):
                    unvisited[neighbor] = new_distance
                    parents[neighbor] = current_node

            visited[current_node] = current_distance
            unvisited.pop(current_node)

            if current_node == end:
                break

        path = []
        node = end
        while parents[node] is not None:
            for link in self.Links:
                if link.O_link.ID == parents[node] and link.D_link.ID == node:
                    path.insert(0, link.ID)
            node = parents[node]

        return path

    def all_none(self):  # 全有全无分配函数
        all_none_linkflow = [0 for i in range(len(self.Links))]
        for i in range(len(self.Links)):
            self.Links[i].Traveltime = self.Links[i].FFT * (
                        1 + self.Links[i].B * (float(self.Linkflows[i]) / float(self.Links[i].Capacity)) ** self.Links[
                    i].Power)  # 更新路段行驶时间
            all_none_linkflow[i] = 0
        for i in range(len(self.Origins)):
            o_node = self.Origins[i].origin_node.ID
            for j in range(len(self.Origins[i].Destination)):
                d_node = self.Origins[i].Destination[j]
                demand = self.Origins[i].odlinks_demand[j]
                Djpathlink = self.Dijkstra_path(o_node, d_node)  # 找最短路
                for index in Djpathlink:
                    all_none_linkflow[index - 1] += demand  # 将流量加载到最短路上
        return all_none_linkflow

    def getUEerr(self):  # 计算UE误差
        sum1 = 0
        for i in range(len(self.Links)):
            self.Links[i].Traveltime = self.Links[i].FFT * (
                    1 + self.Links[i].B * (self.Linkflows[i] / self.Links[i].Capacity) ** self.Links[i].Power)
            sum1 += self.Links[i].Traveltime * self.Linkflows[i]  # 计算流量与行驶时间的乘积（UE公式中的积分项）
        sum2 = 0
        for i in range(len(self.Origins)):
            for j in range(len(self.Origins[i].Destination)):
                demand = self.Origins[i].odlinks_demand[j]
                cost = self.Dijkstra_path(self.Origins[i].origin_node.ID, self.Origins[i].Destination[j])
                sum2 += demand * cost[-1]  # 计算需求与行驶时间的乘积，使用路径上的最后一个节点的成本
        return 1 - sum2 / sum1

    def Optfunction(self, Descent, Lamuda):  # 计算函数值，用于一维搜索
        Sum = 0
        for i in range(len(self.Links)):
            x = self.Linkflows[i] + Lamuda * Descent[i]
            Sum += Descent[i] * (
                        self.Links[i].FFT * (1 + self.Links[i].B * (x / self.Links[i].Capacity) ** self.Links[i].Power))
        return Sum

    def Frank_Wolfe(self):  # F定义rank-Wolfe主函数
        iter = 0  # 迭代次数
        self.Linkflows = [0 for i in range(len(self.Links))]
        self.Linkflows = self.all_none()
        while self.err > self.max_err:
            oldlinkflow = self.Linkflows
            newlinkflow = self.all_none()
            Descent = []
            for i in range(len(self.Links)):
                Descent.append(newlinkflow[i] - self.Linkflows[i])
            Lamuda = 0
            left = 0  # 二分法求最优步长
            right = 1
            mid = 0
            f_left = self.Optfunction(Descent, left)
            f_right = self.Optfunction(Descent, right)
            f_mid = 0
            if f_left * f_right > 0:
                if abs(f_left) > abs(f_right):
                    Lamuda = right
                else:
                    Lamuda = left

            else:
                while right - left > self.max_err:
                    mid = (left + right) / 2
                    f_left = self.Optfunction(Descent, left)
                    f_right = self.Optfunction(Descent, right)
                    f_mid = self.Optfunction(Descent, mid)
                    if f_left * f_mid > 0:
                        left = mid
                    else:
                        right = mid
                Lamuda = (left + right) / 2

            for i in range(len(self.Links)):  # 更新路段流量
                self.Linkflows[i] += Lamuda * Descent[i]
                print(self.Linkflows)
            iter += 1
            self.err = self.getUEerr()
            # print(self.Linkflows)


net = network()
# 读取路网文件
net.read_nodes('节点.csv')
net.read_link('路段.csv')
net.read_od('OD需求.csv')
net.Frank_Wolfe()
print(net.Linkflows)


