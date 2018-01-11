Optics Expert 사용자 메뉴얼
===========================

소개
----

OpticsExpert는 python으로 만든 빔 광학 계산 프로그램이다. 오래전에 만들어진 빔 광학 프로그램들의 단점(Text User Interface)을 개선하여 편의성을 증가시키는것이 이 프로그램을 만든 이유이다. 이 프로그램은 GPL 3.0 라이선스를 따른다. 따라서 교육용, 상업용을 포함한 모든 목적의 사용이 가능하다.

OpticsExpert는 GitHub에서 관리와 업데이트가 공식적으로 이루어진다. 주소는 다음과 같다. http://github.com/garamhahn/OpticsExpert

이 문서는 프로그램의 설치, 이용 방법, 원리를 모두 다룬다. 설치는 설치 섹션에, 원리는 부록에, 이용 방법은 나머지 섹션에 설명하였다.

이 프로그램은 2017년 대한민국 교육부(이공학개인기초연구과제, NRF-2017R1D1A1B03036184)의 지원을 받아 제작되었다.

설치
----

### Python

OpticsExpert를 실행하려면 python이 PC에 설치되어 있어야 합니다. 개발과 테스트는 python 2.x 버전에서 이루어졌습니다. python 3.x 이용자분들께서는 version 2.x로 다운그레이드 해서 이용하시거나 가상 환경에서 이용하시는 것을 추천합니다. Python의 공식 홈페이지는 다음과 같습니다. (http://python.org\)

파이선의 설치는 다음의 문서를 참고한다. http://docs.python-guide.org/en/latest/starting/installation/

### 필수 라이브러리들

OpticsExpert는 다음의 세 가지 라이브러리들을 사용한다. 1. matplotlib - 플롯을 그리고 업데이트 하기 위해 1. scipy - 다변수 해석 등의 수치해석 알고리즘을 이용하기 위해 1. wxpython - 창, 버튼, 멀티프로세싱 등을 사용하기 위해

따라서 python의 설치 이후, `pip`을 이용하여 아래와 같은 명령으로 필수 라이브러리들을 설치한다. 나머지 의존성은 `pip`이 모두 알아서 설치할 것이다.

```
$ pip install -U matplotlib scipy wxpython
```

올바로 설치가 되었는지 확인하기 위해 다음과 같은 명령을 사용할 수 있다. `pip list` 옵션은 현재 설치되어있는 라이브러리를 모두 보여준다.

```
$ pip list
.
.
matplotlib (2.1.1) ............ [확인]
..
numpy (1.13.3)     ............ [확인]
.
<중간중간생략>
.
scipy (1.0.0)      ............ [확인]
.
.
wxPython (4.0.0b2) ............ [확인]
.
Garam@NEBUCHADNEZZAR C:\Users\Garam
$
```

그 밖의 `pip`의 사용법은 다음의 문서를 참고한다. https://pip.pypa.io/en/stable/user_guide/

### 소스코드 복사

OpticsExpert는 github를 이용하여 개발/관리된다. 따라서 최신 버전의 소스코드를 얻기 위해서는 다음과 같이 실행한다.\` $ git clone http://github.com/garamhahn/OpticsExpert 위 명령 정상적으로 종료되면, OpticsExpert 소스코드가 동일한 이름의 디렉토리 안에 복사 되어있을 것이다.

둘러보기
--------

설치 섹션의 두 단계 (1) 파이선 설치, (2) 필수 라이브러리 설치, (3) OpticsExpert 다운로드가 잘 되었다면, 이제 코드를 실행할 준비가 되었다. 코드는 OpticsExpert가 복사된 디렉토리에서 아래와 같이 실행한다.

```
$ python OpticsExpert.py
```

![](FIRST_RUN.PNG)

### 메인 윈도우

#### 구조

#### 단축 아이콘

#### 메뉴 설명

### 입력 창

#### 구조

#### 메뉴 설명

### 출력 창

#### 구조

#### 메뉴 설명

빔 정의
-------

빔라인 디자인
-------------

트랙킹
------

기타 기능들
-----------

부록
----
