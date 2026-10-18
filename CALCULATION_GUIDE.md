# Channel Model Generator — 계산 방식과 사용법

이 프로그램은 초기 설계 및 채널 모델 실험용입니다. 형상 모델은 2D 근사,
via는 등가회로 근사입니다. 실제 PCIe Gen5+/DP/eDP 적합성 판정을 수행하거나
3D EM solver를 대체하지 않습니다. 측정/EM S-parameter 파일을 블록별로
사용하면 해당 데이터의 결합·반사 특성을 그대로 네트워크에 포함합니다.

## 1. 실행 및 파일 구성

```powershell
# 현재 PC: 번들 Python을 찾는 실행 스크립트
& .\run_channel_generator.ps1

# 일반 Python 설치 환경
python -m pip install -r requirements.txt
python Channel_model_generator.py

# GUI 없이 프로젝트 → Touchstone
python Channel_model_generator.py --project examples\line_via_line.json --output output\channel.s8p

# 검증
python Channel_model_generator.py --self-test
```

출력 폴더는 미리 만들거나 기존 폴더를 지정하세요. GUI 저장 대화상자에서
선택한 파일 확장자는 실제 포트 수에 맞게 수정됩니다.

- `Channel_model_generator.py`: GUI, CLI, 프로젝트 저장, 그래프, 결과 출력
- `channel_engine.py`: 물리 선로, 측정 길이 변환, via, Touchstone, 네트워크 합성
- `test_channel_engine.py`: 수치/물리 불변량 테스트
- `model_library/`: 사용자 부품 프리셋 레지스트리와 제조사 자료 안내
- `examples/`: 초기 선로, 여러 lane + via, 측정 기준 파일과 길이 변환 예제

Tkinter는 일반 Windows Python 설치의 Tcl/Tk 옵션이 필요합니다. 현재 Codex
번들 런타임은 샌드박스 안에서 Tcl 파일을 읽지 못했지만, 정상 사용자 실행
환경에서는 숨겨진 창으로 GUI 초기화·블록 편집·계산·그래프를 검사했습니다.

## 2. 사용자 흐름

1. 처음에는 FR4 예시 / differential 1 pair / 50 mm 선로가 생성됩니다.
2. 시작할 때만 선폭을 계산해 평가 주파수에서 Re(Zdiff)=100 Ω에 맞춥니다.
3. Single을 선택하면 현재 형상 선로를 50 Ω로 초기화합니다.
4. 이후 선폭/간격/재질 변경은 형상을 유지하고 상단 계산 임피던스만 갱신합니다.
5. `목표 임피던스로 재설계`를 누를 때만 선폭을 다시 풉니다.
6. 블록 추가/순서 변경 후 계산합니다. 계산 이후 입력 변경 시 재계산 전에는
   이전 결과를 내보낼 수 없습니다.
7. Touchstone, 동일 이름의 `.project.json`, `.report.json`을 함께 저장합니다.

목표값은 선로의 특성 임피던스, `zref_ohm`은 S-parameter power-wave 기준입니다.
두 값은 다릅니다. 포트 reference를 50 Ω로 변경한다고 실제 선로가 50 Ω로
바뀌지 않습니다. 인터페이스 선택은 현재 **프로젝트 태그**입니다. 세대별
권장 임피던스·분석 대역을 자동 인증/강제하지 않습니다.

## 3. 포트 정의

Single-ended 출력 순서:

```
IN_L1_P, IN_L1_N, IN_L2_P, IN_L2_N, ...,
OUT_L1_P, OUT_L1_N, OUT_L2_P, OUT_L2_N, ...
```

Single 신호에서는 P/N 없이 IN_L1,...,OUT_L1,...입니다.
Differential N pair는 4N ports, single N 경로는 2N ports입니다.
여기서 lane은 **한 방향 pair**입니다. PCIe 논리 lane은 TX와 RX 두 pair이므로
필요한 방향을 별도로 모델링하거나 pair 수를 맞춰야 합니다.

`port_order`는 파일의 물리 포트를 내부 순서로 나열한 **1-based permutation**입니다.
예: 원본이 [IN_P, OUT_P, IN_N, OUT_N]이면 `[1,3,2,4]`.
파일 선택 시 `[1,2,...]`를 제안하므로 실제 제조사 매핑을 확인해야 합니다.

## 4. 형상 기반 선로: `capacitance_air`, `line_modes`, `physical_line`

### 4.1 단면과 정전용량

모든 도체가 같은 층에 평행하게 놓인 직사각형 단면을 가정합니다.
선폭 w, 두께 t, 하부 ground와 동박 하단 간 h를 입력합니다.
Stripline은 동박 상단과 상부 ground 간 높이도 지정합니다.
Pair 내부 간격과 lane 간 간격은 모두 **edge-to-edge**입니다.

각 직사각형의 위/아래 면은 각 6개, 양 옆은 각 2개 전하 panel로 나눕니다.
패널 전하와 중심점 전위를 연결하는 행렬 P를 만들고, 도체별 1 V 인가 조건을
풀어 전체 도체의 공기 중 정전용량 행렬 C_air [F/m]을 얻습니다.

- Microstrip ground-plane image kernel:
  `Pij = log(r_image² / r_direct²) / (4π ε0)`
- Stripline parallel-ground Green kernel:
  `Pij = log[(cosh(π dx/H)-cos(π(y_i+y_j)/H)) /
             (cosh(π dx/H)-cos(π(y_i-y_j)/H))] / (4π ε0)`
- Self panel은 평균 log 특이항을 `ds/(2e)`의 등가 거리로 처리합니다.
- 도체-패널 소속 행렬을 B라 하면 `C_air = Bᵀ solve(P,B)`입니다.

이 계산은 중심점 collocation과 self-panel 근사를 사용합니다. 실제 제조
단면에 대해 mesh 수렴/상용 solver 비교를 완료한 field solver는 아닙니다.

### 4.2 재질과 손실

FR4/low-loss/PTFE-like는 제조사 인증 소재가 아닌 수정 가능한 예시값입니다.
Stripline은 균일 유전체, microstrip은 단면 전체에 공통인 유효 유전율을
사용합니다. Microstrip fill factor는 w/h 기반 경험적 근사입니다.

분산 옵션을 켜면 e^(jωt) convention에서 다음 유전체 함수를 사용합니다.

```
q(f) = log((f_high + j f)/(f_low + j f))
strength = Dk_ref * Df_ref / (-Im(q(f_ref)))
eps_inf = Dk_ref - strength * Re(q(f_ref))
eps_r(f) = eps_inf + strength * q(f)
```

이는 유한 주파수 대역의 연속 Debye pole 분포 형태입니다. 기준 주파수에서
입력 Dk와 Df를 맞춥니다. `eps_inf < 1` 같은 부적절한 입력은 거부합니다.
분산을 끄면 상수 Dk/Df이며 광대역 시간영역 인과성을 보장하지 않습니다.

```
L = μ0 ε0 inverse(C_air)
C(f) = Re(eps_eff(f)) C_air
G(f) = -ω Im(eps_eff(f)) C_air
Rdc = 1/(σ w t)
Rskin = sqrt(π f μ0 / σ)/(2(w+t))
roughness_factor = 1 + 2/π atan(1.4 (roughness/skin_depth)^2)
R(f) = sqrt(Rdc² + (Rskin*roughness_factor)²) I
```

도체 skin loss는 크기 기반 보간 근사입니다. 내부 인덕턴스 분산·ground
return loss·proximity effect를 완전히 풀지 않았으므로 **유전체의 인과적
함수 사용이 전체 채널의 인과성을 보증하지는 않습니다.**

### 4.3 전송선 S-parameter

C_air의 고유벡터 U로 modal RLGC를 대각화합니다. 각 mode에 대해:

```
gamma = sqrt((R+jωL)(G+jωC))
Zc = sqrt((R+jωL)/(G+jωC))
p = exp(-gamma * length)
rho = (Zc-Zref)/(Zc+Zref)
S11_mode = rho*(1-p²)/(1-rho² p²)
S21_mode = p*(1-rho²)/(1-rho² p²)
```

고유벡터 U로 도체 좌표계에 되돌려 2N-conductor terminal S를 구성합니다.
길이 0에서는 ideal through입니다. 지수 감쇠 형태를 사용해 긴 선로의
cosh/sinh overflow를 줄였습니다.

### 4.4 임피던스 표시 및 자동 선폭

`Z_matrix = U diag(Zc) Uᵀ`.
Single 표시값은 Z_matrix의 대각, differential pair의 값은
`Zpp + Znn - 2 Zpn`입니다. 다중 도체의 특성 임피던스 **행렬**에 기반한
값이므로 다른 도체의 특정 실제 종단에서 측정한 입력 임피던스와 다릅니다.

선폭 합성은 평가 주파수에서 lane 평균 `Re(Z)`와 목표값 차이를 이분법으로
풉니다. 시작 탐색 범위는 0.005 mm부터 높이에 비례한 최대 폭입니다.
동일 선폭인 여러 lane에서는 바깥 lane과 안쪽 lane의 임피던스가 다를 수
있습니다. 평균을 맞추고 lane별 값을 표시합니다. 개별 lane 선폭 합성은
현재 지원하지 않습니다.

### 4.5 Crosstalk 의미

전하 행렬의 비대각 성분에서 도체 결합을 계산하므로 임의의 crosstalk dB를
입력하지 않습니다. 그러나 microstrip에서 mode마다 다른 유전체 filling과
속도 차이를 풀지 않으므로 **특히 FEXT의 정량 정확도가 제한됩니다.**
단면은 평행하고 균일하며 ground plane은 이상적이라고 가정합니다.
Solder mask, weave, radiation, 층별 비균일 재질은 미지원입니다.

`extra_lengths_mm`는 공통 결합 구간 뒤에 붙는 **독립 tail**입니다. 같은
구간에서 lane별 길이만 임의로 달라지는 geometry를 해석하는 기능은 아닙니다.

## 5. 측정 기준 선로: `resize_uniform`

필수: single-ended S 파일, 포트 매핑, 기준 길이, 새 길이,
`uniform_line_confirmed=true`.

파일은 균일 선로의 양 끝 기준면이어야 합니다. Connector/via/fixture가
남아 있으면 그 부분까지 길이 배율이 적용되므로 사용하면 안 됩니다.
현 버전은 자동 fixture de-embedding 또는 두 길이 측정 추출을 수행하지 않습니다.

양쪽 conductor group으로 분할한 S 행렬에서 wave-chain T를 만듭니다.

```
[a_left; b_left] = T [b_right; a_right]
T11 = inverse(S21)
T12 = -inverse(S21) S22
T21 = S11 inverse(S21)
T22 = S12 - S11 inverse(S21) S22

T_new = exp(log(T_reference) * new_length/reference_length)
```

행렬 함수는 NumPy 고유값 분해로 계산합니다. 원본 주파수축에서 고유벡터
overlap/고유값 거리로 mode를 추적하고 고유값 위상을 unwrap한 후 배율을
적용합니다. 마지막에 출력 주파수축으로 보간합니다. 따라서 단순 S21 dB
배율과 달리 선로 특성 임피던스 불일치에 따른 반사도 길이 변화에 반영됩니다.

검사:

- S21 역행렬과 고유벡터 조건수가 지나치게 나쁘면 거부
- 첫 고유위상이 큰 경우 낮은 시작 주파수 요청
- 큰 주파수별 위상 변화 거부
- 기준 길이 재생성 오차 검사
- 길이 변환 결과의 수동성 검사

**한계:** 첫 주파수에서 이미 여러 회전한 위상은 작은 위상으로 보일 수
있습니다. 데이터만으로 branch ambiguity를 완전히 탐지할 수 없습니다.
기준 선로가 첫 주파수에서 충분히 짧고, 대략 `Δf * delay << 0.5`인 조밀한
샘플링인지 사용자가 확인해야 합니다. 재생성 오차가 작다고 물리적으로
올바른 균일 선로 모델이라는 뜻은 아닙니다. 축퇴/교차 mode가 있는 다중
lane 길이 변환은 실험적이며 독립 검증이 필요합니다.

새 길이가 기준 길이와 같으면 원본을 사용합니다. 새 길이 0은 through입니다.
파일 모드에서 원래 선폭/Dk를 역추정해 편집하거나 자동으로 100 Ω에 맞추지 않습니다.

## 6. Connector / 측정 Via: `imported_block`

파일은 고정 블록으로 사용합니다. 포트 순서를 바꾼 다음 공통 reference로
renormalize하고 공통 주파수축에 보간합니다. 파일 안의 다중 lane 결합과
mode conversion을 유지합니다.

```
r = (Znew-Zold)/(Znew+Zold)
Snew = (S-rI) inverse(I-rS)
```

보간은 real/imag 선형이며 파일 범위 밖으로 외삽하지 않습니다. 너무 성긴
데이터의 긴 지연 위상을 복구하지 못하므로 원본의 충분한 샘플링이 필요합니다.

1-pair 파일 복제는 사용자가 명시적으로 켤 때만 허용합니다. 이 경우 lane
간 결합은 **모델에 없는 것**이며 물리적인 crosstalk=0을 검증한 것이 아닙니다.
출력 notes에 이 제한이 남습니다.

현재 full-matrix Touchstone 1.x S/RI, MA, DB, 동일한 실수 reference만 지원합니다.
Touchstone 2.0, unequal/complex references, mixed-mode 입력, noise section은
추측해서 읽지 않고 거부합니다. Export 시 single-ended 1.x로 바꿔 넣으세요.
Two-port는 S11,S21,S12,S22, 3-port 이상은 row-major 표준 순서를 사용합니다.

모델의 기준면이 다른 블록과 맞는지는 includes 설명을 통해 사람이 확인해야
합니다. 모델 설명만으로 중복 pad/via를 자동 탐지하는 기능은 없습니다.
사용하지 않는 포트를 임의로 삭제하거나 자동 종단하지 않습니다. 필요한
종단이 반영된 reduced network를 입력해야 합니다.

## 7. Via 근사: `via_network`

Through/backdrilled/blind/microvia 선택은 **치수 초기값**만 바꿉니다.
모두 동일한 근사 등가회로 엔진을 사용합니다. 제조 공정 인증값이 아닙니다.

- 사용 barrel: 직렬 R + jωL
- Pad/antipad: shunt C
- 남은 stub: open-ended transmission-line admittance
- Pair 상호 인덕턴스: 사용자 가정 k*L

길이/직경을 inch로 환산한 경험식:

```
L[nH] = 5.08 h_in (ln(4h/d)+1)
C[pF] = 1.41 Dk h_in pad_in/(antipad_in-pad_in)
Y_stub = tanh(gamma_stub * stub_length)/Z_stub
```

전체 shunt C/stub admittance의 절반을 직렬 RL 양쪽에 배치한 대칭 π 근사입니다.
Stub propagation은 sqrt(Dk)/c 속도와 작은 감쇠 0.001 Np/m를 사용합니다.
Stub Z는 형상에서 엄밀히 추출하지 않고 사용자 값입니다.

**Ground via 배치, plane cavity, 인접 via array crosstalk, 3D mode 변환은
계산하지 않습니다.** 현재 stackup layer table이 없으므로 사용 길이/stub를
자동 층 계산하지 않습니다. 이런 효과가 중요한 Gen5+ 채널은 via 영역 전체의
EM/측정 S-parameter로 대체하세요. 짧은 via도 모델 검증 없이 정확하다고
가정하지 않습니다.

## 8. 전체 채널 합성: `cascade`

앞 블록 A의 오른쪽과 뒷 블록 B의 왼쪽을 연결하는 Redheffer 연결을 사용합니다.

```
X = solve(I-A22 B11, A21)
Y = solve(I-A22 B11, A22 B12)
S11 = A11 + A12 B11 X
S12 = A12 (B12+B11 Y)
S21 = B21 X
S22 = B22+B21 Y
```

따라서 블록 간 반사와 기존 다중 포트 결합을 유지합니다. 블록 A와 블록 B
사이의 포트 외 **공간적인 근접 결합**을 새로 계산하지는 않습니다. 가까운
connector/launch/via가 서로 영향을 주면 하나의 결합 EM 모델로 묶어야 합니다.

각 주파수에서 최대 singular value <=1인지 검사합니다. 최종 계산은 1.005를
넘으면 거부하며, 허용된 수치 오차와 실제 최대값은 report에 기록합니다.
수동성 강제 보정과 완전한 인과성 인증은 제공하지 않습니다.

## 9. 결과 그래프

Differential에서는 각 side의 P/N power wave를 직교 변환합니다.

```
ad = (ap-an)/sqrt(2)
ac = (ap+an)/sqrt(2)
```

Differential physical reference는 2*Zref, common은 Zref/2입니다.
그래프는 변환된 Sdd/SCD를 사용하지만 파일 출력은 single-ended입니다.

- Insertion loss: 20log10|Sdd21| (음수 dB convention)
- Return loss: 20log10|Sdd11| (양의 return-loss 값이 아닌 S11 dB)
- NEXT/FEXT: lane 1을 aggressor로 다른 lane으로 전달되는 차동 성분
- D to C: 각 lane differential 입력 → 동일 lane common 출력
- Group delay: -d unwrap(arg(S21))/dω, ps

그래프 바닥은 -160 dB입니다. Transmission null 근처 group delay는 불안정할
수 있습니다. 전체 채널 또는 개별 블록을 선택해 비교할 수 있습니다.
Eye/BER/규격 mask, TDR, DC 외삽, mixed-mode 파일 export는 현재 미구현입니다.

## 10. 실제 검증 범위

자동 테스트는 초기 50/100 Ω 합성, 선폭-임피던스 방향성, 길이 0,
두 구간 합성과 전체 길이 일치, 물리 선로로 생성한 기준 파일 길이 변환,
다중 lane 길이 변환/결합/상호성/수동성, 재질/stripline, Zref 왕복 변환,
Touchstone 순서·포맷·왕복, 사용자 포트 재배열, 독립 복제,
독립 추가 길이, via 포함 채널 수동성, 잘못된 입력 차단을 검사합니다.

이 테스트들은 엔진의 수치 일관성 검증입니다. 실측 보드·상용 solver와의
정량 오차 검증, 특정 제조사 connector 검증, PCIe/DP 인증 시험은 아닙니다.
