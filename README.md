# Channel Model Generator

`Channel_model_generator.py`를 실행하는 Tkinter GUI + NumPy 채널 계산 도구입니다.

```powershell
& .\run_channel_generator.ps1
```

일반 Python 환경에서는 `pip install -r requirements.txt` 후
`python Channel_model_generator.py`로 실행하세요. Tkinter가 필요합니다.

## 구현된 기능

- 형상 기반 single/differential 선로, 50/100 Ω 초기 선폭 합성
- 선폭/간격 수정 시 형상을 유지하며 임피던스 갱신
- 여러 lane의 2D 결합 근사, 공통 결합 길이와 lane별 독립 tail
- 동일 단면 기준 S-parameter의 길이 변환
- Connector/via S-parameter 블록 삽입, 순서 변경, 파일 프리셋 등록
- Via의 수정 가능한 RL/C/stub 근사 프리셋
- 전체 다중 포트 합성, 삽입/반사/NEXT/FEXT/mode conversion/group delay 그래프
- Single-ended Touchstone 출력, 프로젝트/계산 가정 저장
- GUI 없는 CLI와 자체 테스트

## 먼저 열어볼 예제

- `examples/default_line.json`: 자동 100 Ω 초기 선로
- `examples/line_via_line.json`: 2 differential pair, 선로–via–선로
- `examples/measured_resize_example.json`: 50 mm 기준 파일 → 80 mm

`synthetic_reference_50mm.s4p`는 엔진으로 만든 **합성 검증 데이터**입니다.
실제 보드 실측 데이터가 아닙니다. 생성된 예제 `.s8p`도 실제 제품 검증값이 아닙니다.

Connector 제조사 파일은 현재 번들하지 않았습니다. 파일 블록을 설정한 후
`현재 파일 블록을 프리셋 저장`으로 등록하면 다시 선택해서 사용할 수 있습니다.

## 계산 범위

계산 공식, 코드 함수별 설명, 단위, 포트 매핑, 실측 길이 변환 조건,
지원하지 않는 항목은 **[CALCULATION_GUIDE.md](CALCULATION_GUIDE.md)**에 있습니다.

특히 형상 기반 microstrip crosstalk와 via는 초기 설계용 근사입니다.
PCIe Gen5+/DP/eDP 규격 인증이나 3D EM 정확도를 주장하지 않습니다.
Touchstone 1.x full-matrix single-ended 입력을 지원하며 2.0/mixed-mode 입력은
사전 변환이 필요합니다. 인터페이스 선택은 태그이며 규격 mask는 미구현입니다.
