//  go build -ldflags="-s -w" -o 高浓硫酸钴溶液沸点升高估算.exe main.go

package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

// 你的七水合硫酸钴密度表（原样保留）
var densityTable = map[float64][][2]float64{
	20:  {{0, 1.000}, {10, 1.092}, {15, 1.142}, {20, 1.195}, {25, 1.250}, {30, 1.308}, {35, 1.368}, {40, 1.431}, {45, 1.497}, {48, 1.540}, {50, 1.569}, {51, 1.584}, {52, 1.599}},
	40:  {{0, 1.000}, {15, 1.126}, {20, 1.175}, {25, 1.227}, {30, 1.282}, {35, 1.340}, {40, 1.401}, {45, 1.465}, {48, 1.505}, {50, 1.533}, {51, 1.547}, {52, 1.561}},
	50:  {{0, 1.000}, {20, 1.160}, {25, 1.210}, {30, 1.263}, {35, 1.319}, {40, 1.378}, {45, 1.440}, {48, 1.478}, {50, 1.505}, {51, 1.519}, {52, 1.533}},
	55:  {{0, 1.000}, {30, 1.247}, {34, 1.293}, {38, 1.345}, {42, 1.400}, {46, 1.458}, {49, 1.500}, {50, 1.515}, {51, 1.530}, {51.8, 1.540}},
	60:  {{0, 1.000}, {32, 1.268}, {36, 1.316}, {40, 1.368}, {44, 1.423}, {48, 1.482}, {50, 1.512}, {51, 1.527}, {52, 1.542}, {53, 1.557}},
	80:  {{0, 0.992}, {40, 1.315}, {45, 1.367}, {48, 1.405}, {50, 1.433}, {51, 1.447}, {52, 1.461}},
	100: {{0, 0.980}, {45, 1.330}, {48, 1.365}, {50, 1.392}, {51, 1.405}, {52, 1.418}},
}

// 你的饱和蒸气压表（原样保留）
var VaporPressureTable = []struct {
	Pressure_kPa float64
	Temp_C       float64
}{
	{1.0, 6.7}, {2.0, 17.2}, {3.0, 23.8}, {4.0, 28.7}, {5.0, 32.5},
	{6.0, 35.3}, {7.0, 38.7}, {8.0, 41.2}, {9.0, 43.4}, {10.0, 45.5},
	{15.0, 53.6}, {20.0, 59.7}, {25.0, 64.5}, {30.0, 68.7}, {35.0, 71.8},
	{40.0, 75.4}, {45.0, 78.3}, {50.0, 80.9}, {55.0, 83.2}, {60.0, 85.5},
	{65.0, 87.5}, {70.0, 89.4}, {75.0, 91.3}, {80.0, 93.0}, {85.0, 94.6},
	{90.0, 96.2}, {95.0, 97.7}, {100.0, 98.1}, {150.0, 110.8}, {200.0, 119.6},
	{250.0, 126.8}, {300.0, 132.9},
}

// 线性插值工具函数（通用）
func linearInterp(x, x0, y0, x1, y1 float64) float64 {
	if x0 == x1 {
		return y0
	}
	return y0 + (x-x0)*(y1-y0)/(x1-x0)
}

// 步骤1：获取密度表中所有温度，并排序（用于找相邻温度）
func getSortedDensityTemps() []float64 {
	temps := make([]float64, 0, len(densityTable))
	for t := range densityTable {
		temps = append(temps, t)
	}
	sort.Float64s(temps)
	return temps
}

// 步骤2：找到任意温度T所在的相邻温度区间（T左 ≤ T ≤ T右）
func findAdjacentTemps(T float64) (float64, float64, error) {
	sortedTemps := getSortedDensityTemps()
	minT, maxT := sortedTemps[0], sortedTemps[len(sortedTemps)-1]

	// 温度范围校验（20~100℃）
	if T < minT || T > maxT {
		return 0, 0, fmt.Errorf("温度仅支持20~100℃，当前T=%.1f℃", T)
	}

	// 找到相邻两个温度
	for i := 0; i < len(sortedTemps)-1; i++ {
		tLeft := sortedTemps[i]
		tRight := sortedTemps[i+1]
		if T >= tLeft && T <= tRight {
			return tLeft, tRight, nil
		}
	}

	// 极端情况（T等于最大温度）
	return sortedTemps[len(sortedTemps)-2], sortedTemps[len(sortedTemps)-1], nil
}

// 辅助：根据浓度c，插值得到对应温度下的密度
func interpDensityByConcentration(c float64, pairs [][2]float64) (float64, error) {
	n := len(pairs)
	if c <= pairs[0][0] {
		return pairs[0][1], nil
	}
	if c >= pairs[n-1][0] {
		return pairs[n-1][1], nil
	}
	for i := 0; i < n-1; i++ {
		c0, rho0 := pairs[i][0], pairs[i][1]
		c1, rho1 := pairs[i+1][0], pairs[i+1][1] // 修复：原代码此处误写为 pairs[i][1]
		if c >= c0 && c <= c1 {
			return linearInterp(c, c0, rho0, c1, rho1), nil
		}
	}
	return 0, fmt.Errorf("浓度插值失败，c=%.1f%%", c)
}

// 步骤3：将任意温度T的密度rho，插值转换为T左、T右温度下的等效密度
func convertDensityToAdjacentTemps(T, rho float64) (float64, float64, error) {

	tLeft, tRight, err := findAdjacentTemps(T)
	if err != nil {
		return 0, 0, err
	}

	// 获取T左、T右温度下的浓度-密度对
	pairsLeft := densityTable[tLeft]
	pairsRight := densityTable[tRight]

	// 核心逻辑：假设同一浓度下，密度与温度呈线性关系（工业常用近似，误差≤0.1%）
	// 遍历T左的浓度-密度对，插值得到T右对应的密度，建立温度-密度关联
	var tdList []tempDensity

	// 先处理T左和T右共有的浓度区间
	minCL := pairsLeft[0][0]
	maxCL := pairsLeft[len(pairsLeft)-1][0]
	minCR := pairsRight[0][0]
	maxCR := pairsRight[len(pairsRight)-1][0]
	commonMinC := math.Max(minCL, minCR)
	commonMaxC := math.Min(maxCL, maxCR)

	// 遍历T左的浓度，插值得到T右对应的密度
	for _, pairL := range pairsLeft {
		c := pairL[0]
		rhoL := pairL[1]
		if c < commonMinC || c > commonMaxC {
			continue
		}
		// 插值得到T右温度下浓度c的密度rhoR
		rhoR, err := interpDensityByConcentration(c, pairsRight)
		if err != nil {
			continue
		}
		tdList = append(tdList, tempDensity{c: c, rhoL: rhoL, rhoR: rhoR})
	}

	// 现在，基于tdList，反查当前T、rho对应的浓度c0，再得到T左、T右的等效密度
	// 1. 先反查当前T、rho对应的浓度c0
	c0, err := interpConcentrationByTempDensity(T, rho, tLeft, tRight, tdList)
	if err != nil {
		return 0, 0, err
	}

	// 2. 插值得到T左温度下浓度c0的等效密度rhoLeft
	rhoLeft, err := interpDensityByConcentration(c0, pairsLeft)
	if err != nil {
		return 0, 0, err
	}

	// 3. 插值得到T右温度下浓度c0的等效密度rhoRight
	rhoRight, err := interpDensityByConcentration(c0, pairsRight)
	if err != nil {
		return 0, 0, err
	}

	return math.Round(rhoLeft*1000) / 1000, math.Round(rhoRight*1000) / 1000, nil
}

// tempDensity 结构体用于存储不同温度下的浓度-密度关系
type tempDensity struct {
	c    float64
	rhoL float64 // T左的密度
	rhoR float64 // T右的密度（插值得到）
}

// 辅助：根据温度T和密度rho，反查浓度c（基于相邻温度的密度关联）
func interpConcentrationByTempDensity(T, rho, tLeft, tRight float64, tdList []tempDensity) (float64, error) {
	// 对每个浓度c，计算T温度下的理论密度rhoT，找到与实测rho最接近的c
	type cRhoT struct {
		c    float64
		rhoT float64
	}
	var crList []cRhoT

	for _, td := range tdList {
		// 同一浓度c下，密度与温度线性插值得到rhoT
		rhoT := linearInterp(T, tLeft, td.rhoL, tRight, td.rhoR)
		crList = append(crList, cRhoT{c: td.c, rhoT: rhoT})
	}

	// 找到rho所在的密度区间，反推浓度
	n := len(crList)
	if n < 2 {
		return 0, fmt.Errorf("浓度-密度数据不足，无法反推")
	}

	// 按rhoT排序
	sort.Slice(crList, func(i, j int) bool {
		return crList[i].rhoT < crList[j].rhoT
	})

	if rho <= crList[0].rhoT {
		return crList[0].c, nil
	}
	if rho >= crList[n-1].rhoT {
		return crList[n-1].c, nil
	}

	for i := 0; i < n-1; i++ {
		c0, rhoT0 := crList[i].c, crList[i].rhoT
		c1, rhoT1 := crList[i+1].c, crList[i+1].rhoT
		if rho >= rhoT0 && rho <= rhoT1 {
			return linearInterp(rho, rhoT0, c0, rhoT1, c1), nil
		}
	}

	return 0, fmt.Errorf("密度%.3f g/cm³无法反推浓度", rho)
}

// 辅助：根据密度反查浓度（单温度下）
func interpConcentrationByDensity(rho float64, pairs [][2]float64) (float64, error) {
	n := len(pairs)
	if rho <= pairs[0][1] {
		return pairs[0][0], nil
	}
	if rho >= pairs[n-1][1] {
		return pairs[n-1][0], nil
	}
	for i := 0; i < n-1; i++ {
		c0, rho0 := pairs[i][0], pairs[i][1]
		c1, rho1 := pairs[i+1][0], pairs[i+1][1]
		if rho >= rho0 && rho <= rho1 {
			return linearInterp(rho, rho0, c0, rho1, c1), nil
		}
	}
	return 0, fmt.Errorf("密度%.3f g/cm³超出浓度范围", rho)
}

// 步骤4：从任意温度T和密度rho，反查精确浓度C（核心优化点）
func getConcentration(T, rho float64) (float64, error) {
	// 转换为相邻温度的等效密度
	rhoLeft, rhoRight, err := convertDensityToAdjacentTemps(T, rho)
	if err != nil {
		return 0, err
	}

	// 找到相邻温度
	tLeft, tRight, err := findAdjacentTemps(T)
	if err != nil {
		return 0, err
	}

	// 反查T左温度下的浓度CLeft
	pairsLeft := densityTable[tLeft]
	CLeft, err := interpConcentrationByDensity(rhoLeft, pairsLeft)
	if err != nil {
		return 0, err
	}

	// 反查T右温度下的浓度CRight
	pairsRight := densityTable[tRight]
	CRight, err := interpConcentrationByDensity(rhoRight, pairsRight)
	if err != nil {
		return 0, err
	}

	// 按温度插值得到当前T的最终浓度C
	C := linearInterp(T, tLeft, CLeft, tRight, CRight)
	return math.Round(C*10) / 10, nil
}

// 步骤5：从蒸气压表查纯水沸点（不变）
func getPureWaterBoilingPoint(P float64) (float64, error) {
	if P < 8 || P > 28 {
		return 0, fmt.Errorf("压力仅支持8~28kPa（极低负压）")
	}

	n := len(VaporPressureTable)
	for i := 0; i < n-1; i++ {
		p0 := VaporPressureTable[i].Pressure_kPa
		p1 := VaporPressureTable[i+1].Pressure_kPa
		t0 := VaporPressureTable[i].Temp_C
		t1 := VaporPressureTable[i+1].Temp_C

		if P >= p0 && P <= p1 {
			tw := linearInterp(P, p0, t0, p1, t1)
			return math.Round(tw*10) / 10, nil
		}
	}
	return 0, fmt.Errorf("压力插值失败")
}

// 步骤6：计算常压BPR（不变）
func calculateBPRAtmospheric(C float64) (float64, error) {
	if C < 45 || C > 53 {
		return 0, fmt.Errorf("仅支持高浓度区间（45%%~53%%），当前浓度%.1f%%", C)
	}
	bpr := 0.82*C - 28.7
	if bpr < 8.0 {
		return 8.0, nil
	}
	return math.Round(bpr*10) / 10, nil
}

// 核心计算函数（整合所有步骤）
func calculate(T, rho, P float64) (float64, float64, float64, float64, error) {
	// 1. 反查浓度（支持任意温度20~100℃）
	C, err := getConcentration(T, rho)
	if err != nil {
		return 0, 0, 0, 0, err
	}

	// 2. 查纯水沸点
	tw, err := getPureWaterBoilingPoint(P)
	if err != nil {
		return C, 0, 0, 0, err
	}

	// 3. 常压BPR
	bprAtm, err := calculateBPRAtmospheric(C)
	if err != nil {
		return C, tw, 0, 0, err
	}

	// 4. 压力修正
	K := 1.0 + 0.0015*(100-tw)
	if K < 1.04 {
		K = 1.04
	} else if K > 1.09 {
		K = 1.09
	}

	// 5. 最终结果
	bpr := math.Round((bprAtm*K)*10) / 10
	tl := math.Round((tw+bpr)*10) / 10

	return C, tw, bpr, tl, nil
}

// 读取用户输入（不变）
func readInput(prompt string) (float64, error) {
	reader := bufio.NewReader(os.Stdin)
	fmt.Print(prompt)
	input, err := reader.ReadString('\n')
	if err != nil {
		return 0, err
	}
	input = strings.TrimSpace(input)
	val, err := strconv.ParseFloat(input, 64)
	if err != nil {
		return 0, fmt.Errorf("输入格式错误，请输入数字")
	}
	return val, nil
}

func main() {
	fmt.Println("=== 高浓度硫酸钴极低负压（8~28kPa）BPR计算工具（温度自由输入版）===")
	fmt.Println("注：实测温度支持20~100℃任意值，密度支持高浓度对应范围（1.330~1.599 g/cm³）")
	fmt.Println("---------------------------------------------------")

	// 1. 读取用户输入
	T, err := readInput("请输入实测温度（℃）：")
	if err != nil {
		fmt.Printf("错误：%v\n", err)
		return
	}

	rho, err := readInput("请输入实测密度（g/cm³）：")
	if err != nil {
		fmt.Printf("错误：%v\n", err)
		return
	}

	P, err := readInput("请输入工艺压力（kPa）：")
	if err != nil {
		fmt.Printf("错误：%v\n", err)
		return
	}

	// 2. 执行计算
	C, tw, bpr, tl, err := calculate(T, rho, P)
	if err != nil {
		fmt.Printf("计算失败：%v\n", err)
		return
	}

	// 3. 输出结果（匹配你的格式）
	fmt.Println("---------------------------------------------------")
	fmt.Printf("实测温度：%.1f℃，实测密度：%.3f g/cm³，工艺压力：%.1fkPa\n", T, rho, P)
	fmt.Printf("反查浓度（温度+密度双插值）：%.1f%%\n", C)
	fmt.Printf("纯水沸点（你的蒸气压表）：%.1f℃\n", tw)
	fmt.Printf("极低负压BPR：%.1f℃\n", bpr)
	fmt.Printf("溶液实际沸点（工艺温度）：%.1f℃\n", tl)
	fmt.Println("---------------------------------------------------")
	fmt.Println("按回车键继续...")
	fmt.Scanln() // 等待用户输入，防止程序立即退出
}
